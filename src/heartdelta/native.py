"""Native-grid longitudinal AHA measurement."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import SimpleITK as sitk
from .aha import LEVELS, NAMES
from .metrics import eqd2_factor, geud


def _same_grid(a: sitk.Image, b: sitk.Image, tolerance: float = 1e-5) -> bool:
    return a.GetSize() == b.GetSize() and all(np.allclose(x, y, atol=tolerance, rtol=0) for x, y in ((a.GetSpacing(), b.GetSpacing()), (a.GetOrigin(), b.GetOrigin()), (a.GetDirection(), b.GetDirection())))


def _labels_on_ct(path: str | Path, ct: sitk.Image) -> np.ndarray:
    labels = sitk.ReadImage(str(path))
    if not _same_grid(labels, ct):
        labels = sitk.Resample(labels, ct, sitk.Transform(), sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
    return sitk.GetArrayFromImage(labels)


def _dose_on_ct(path: str | Path, ct: sitk.Image) -> np.ndarray:
    dose = sitk.ReadImage(str(path))
    if not _same_grid(dose, ct):
        dose = sitk.Resample(dose, ct, sitk.Transform(), sitk.sitkLinear, 0, sitk.sitkFloat32)
    return sitk.GetArrayFromImage(dose).astype(float)


def extract_case(case: dict, alpha_beta: float = 3.0) -> pd.DataFrame:
    """Measure RT dose and native HU at every available timepoint."""
    cid = str(case["case_id"])
    rows: list[dict] = []
    dose_by_segment: dict[int, dict] = {}
    rt_ct_path = case.get("ct_rt") or case.get("baseline_ct")
    rt_aha_path = case.get("aha17_rt") or case.get("aha17")
    if rt_ct_path and rt_aha_path and case.get("dose"):
        ct = sitk.ReadImage(str(rt_ct_path)); labels = _labels_on_ct(rt_aha_path, ct); dose = _dose_on_ct(case["dose"], ct)
        total = float(case.get("total_dose_gy", case.get("total_dose_prior", 50)))
        fractions = float(case.get("fractions", case.get("num_fraction_planned", total / 2)))
        factor = eqd2_factor(total, fractions, alpha_beta)
        spacing = np.prod(ct.GetSpacing()) / 1000.0
        for segment in range(1, 18):
            values = dose[labels == segment]
            dose_by_segment[segment] = {
                "volume_cc": float(len(values) * spacing), "mean_dose_gy": float(np.mean(values)) if len(values) else np.nan,
                "median_dose_gy": float(np.median(values)) if len(values) else np.nan,
                "maximum_dose_gy": float(np.max(values)) if len(values) else np.nan,
                "mean_dose_eqd2": float(np.mean(values) * factor) if len(values) else np.nan,
                **{f"geud_eqd2_a{a}": geud(values * factor, a) for a in (1, 2, 3, 4, 5)},
            }
    for tp in ("rt", "fu1", "fu2", "fu3"):
        ct_path = case.get(f"ct_{tp}") if tp != "rt" else rt_ct_path
        aha_path = case.get(f"aha17_{tp}") if tp != "rt" else rt_aha_path
        if not ct_path or not aha_path or not Path(ct_path).is_file() or not Path(aha_path).is_file():
            continue
        ct = sitk.ReadImage(str(ct_path)); hu = sitk.GetArrayFromImage(ct).astype(float); labels = _labels_on_ct(aha_path, ct)
        for segment in range(1, 18):
            values = hu[labels == segment]
            row = {"case_id": cid, "timepoint": tp, "segment": segment, "segment_name": NAMES[segment], "level": LEVELS[segment], "n_voxels": int(len(values)), "hu_mean": float(np.mean(values)) if len(values) else np.nan, "hu_median": float(np.median(values)) if len(values) else np.nan, "hu_std": float(np.std(values)) if len(values) else np.nan}
            row.update(dose_by_segment.get(segment, {}))
            for field in ("group", "modality", "side", "breathing", "delivery_technique", "target"):
                if field in case: row[field] = case[field]
            rows.append(row)
    frame = pd.DataFrame(rows)
    if frame.empty: return frame
    baseline = frame[frame.timepoint == "rt"][["segment", "hu_mean", "hu_median"]].rename(columns={"hu_mean": "baseline_hu_mean", "hu_median": "baseline_hu_median"})
    frame = frame.merge(baseline, on="segment", how="left")
    frame["delta_hu_mean"] = frame.hu_mean - frame.baseline_hu_mean
    frame["delta_hu_median"] = frame.hu_median - frame.baseline_hu_median
    return frame


def aggregate_levels(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    numeric = [c for c in frame.columns if c not in {"case_id", "timepoint", "segment", "segment_name", "level", "group", "modality", "side", "breathing", "delivery_technique", "target", "n_voxels"} and pd.api.types.is_numeric_dtype(frame[c])]
    for keys, group in frame.groupby(["case_id", "timepoint", "level"], dropna=False):
        row = dict(zip(("case_id", "timepoint", "level"), keys)); row["n_voxels"] = int(group.n_voxels.sum())
        for field in ("group", "modality", "side", "breathing", "delivery_technique", "target"):
            if field in group: row[field] = group[field].iloc[0]
        for col in numeric:
            valid = group[col].notna() & (group.n_voxels > 0)
            row[col] = float(np.average(group.loc[valid, col], weights=group.loc[valid, "n_voxels"])) if valid.any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)
