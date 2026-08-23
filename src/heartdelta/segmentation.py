"""Cardiac chamber and PlatiPy AHA-17 segmentation adapters."""
from __future__ import annotations

import re
from pathlib import Path
import numpy as np
import SimpleITK as sitk

STRUCTURES = ("Atrium_L", "Heart", "Ventricle_L", "Atrium_R", "Ventricle_R", "A_Aorta", "A_LAD")


def _platipy():
    try:
        from platipy.imaging.projects.cardiac.run import run_hybrid_segmentation
        from platipy.imaging.utils.ventricle import generate_left_ventricle_segments
    except ImportError as exc:
        raise RuntimeError("Install the segmentation extra: pip install -e '.[segmentation]'") from exc
    return run_hybrid_segmentation, generate_left_ventricle_segments


def combine_segments(parts: dict[str, sitk.Image]) -> sitk.Image:
    if not parts:
        raise ValueError("AHA generator returned no segments")
    reference = next(iter(parts.values()))
    labels = np.zeros(sitk.GetArrayFromImage(reference).shape, np.uint8)
    for name, image in parts.items():
        match = re.search(r"segment\D*(\d+)", name.lower())
        if match and 1 <= int(match.group(1)) <= 17:
            labels[sitk.GetArrayFromImage(image) > 0] = int(match.group(1))
    result = sitk.GetImageFromArray(labels)
    result.CopyInformation(reference)
    return result


def segment_scan(ct_path: str | Path, output_dir: str | Path, overwrite: bool = False) -> dict[str, str]:
    """Create chamber masks and an AHA-17 map for one native CT."""
    output = Path(output_dir)
    aha_path = output / "aha17.nii.gz"
    if aha_path.is_file() and not overwrite:
        return {"aha17": str(aha_path), **{s: str(output / f"{s}.nii.gz") for s in STRUCTURES if (output / f"{s}.nii.gz").is_file()}}
    output.mkdir(parents=True, exist_ok=True)
    run_hybrid, generate_aha = _platipy()
    ct = sitk.ReadImage(str(ct_path))
    masks, *_ = run_hybrid(ct)
    for name, image in masks.items():
        sitk.WriteImage(image, str(output / f"{name}.nii.gz"))
    required = {name: masks[name] for name in ("Ventricle_L", "Atrium_L", "Ventricle_R", "Heart")}
    aha = combine_segments(generate_aha(required))
    sitk.WriteImage(aha, str(aha_path))
    return {"aha17": str(aha_path), **{s: str(output / f"{s}.nii.gz") for s in STRUCTURES if s in masks}}


def qc_labels(path: str | Path, minimum_voxels: int = 20) -> dict[str, object]:
    labels = sitk.GetArrayFromImage(sitk.ReadImage(str(path)))
    counts = {i: int(np.count_nonzero(labels == i)) for i in range(1, 18)}
    return {
        "labels_present": len([n for n in counts.values() if n]),
        "complete": all(n >= minimum_voxels for n in counts.values()),
        "minimum_segment_voxels": min(counts.values()),
        **{f"segment_{i}_voxels": counts[i] for i in counts},
    }
