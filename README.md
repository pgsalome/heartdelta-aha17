# HeartDelta-AHA17

[![Python](https://img.shields.io/badge/Python-%E2%89%A53.10-3776AB)](pyproject.toml)
[![Tests](https://github.com/pgsalome/heartdelta-aha17/actions/workflows/tests.yml/badge.svg)](https://github.com/pgsalome/heartdelta-aha17/actions/workflows/tests.yml)
[![License](https://img.shields.io/badge/use-noncommercial_only-4B5563)](LICENSE)

HeartDelta-AHA17 generates cardiac structure masks and AHA 17-segment
left-ventricular maps from CT images. It also measures regional radiotherapy
dose and longitudinal CT attenuation changes.

Provide a CT to run segmentation, or a case registry to analyze planning and
follow-up scans. Each CT is measured on its native grid using its own AHA map.

## Table Of Contents

1. [Quick Start](#quick-start)
2. [Segmentation Model](#segmentation-model)
3. [Your Images](#your-images)
4. [Running The Workflow](#running-the-workflow)
5. [Parameters And Outputs](#parameters-and-outputs)
6. [Quality Control](#quality-control)
7. [Tests](#tests)
8. [Citation](#citation)
9. [License](#license)

## Quick Start

Requires Python 3.10 or newer. Linux and an NVIDIA GPU with CUDA-compatible
PyTorch are recommended for automatic segmentation; extraction from reviewed
masks can run on CPU.

To install HeartDelta-AHA17 with the segmentation dependencies:

```bash
git clone --branch main --single-branch https://github.com/pgsalome/heartdelta-aha17.git
cd heartdelta-aha17
python -m venv .venv
source .venv/bin/activate
pip install -e '.[segmentation]'
```

To segment one CT, replace `ct.nii.gz` with your image path and run:

```bash
python - <<'PY'
from heartdelta.segmentation import segment_scan

products = segment_scan("ct.nii.gz", "processed/scan001")
print(products["aha17"])
PY
```

The first run downloads the missing model and atlas files automatically.
Results are saved in `processed/scan001`: `aha17.nii.gz`, `Heart.nii.gz`,
`Ventricle_L.nii.gz`, and the other generated cardiac structure masks.

To regenerate an existing segmentation, pass `overwrite=True` to `segment_scan`.
For analysis with existing reviewed maps, install the base package with
`pip install -e .` instead.

## Segmentation Model

Segmentation uses [PlatiPy's hybrid cardiac workflow](https://github.com/pyplati/platipy/blob/master/platipy/imaging/projects/cardiac/README.md):
a pretrained nnU-Net whole-heart model, cardiac atlas mapping, and geometric
definitions of smaller structures. AHA-17 regions are generated from the
cardiac chamber masks.

No training step or separate HeartDelta checkpoint is required. PlatiPy
downloads the [whole-heart model](https://zenodo.org/record/6585664/files/Task400_OPEN_HEART_3d_lowres.zip?download=1)
and [cardiac atlas](https://zenodo.org/record/6592437/files/open_atlas.zip?download=1)
when they are missing and reuses them on later runs.

| Resource | Default location | Override |
| --- | --- | --- |
| nnU-Net model (`Task400_OPEN_HEART_1FOLD`) | `~/.platipy/nnUNet_models` | `RESULTS_FOLDER` |
| Cardiac atlas | `~/.platipy/cardiac/test_atlas` | `ATLAS_PATH` |

Set the override environment variables before starting Python. The first run
needs internet access to download these resources. For CUDA installation and
segmentation troubleshooting, follow the PlatiPy documentation linked above.

<details>
<summary>Install downloaded model and atlas archives</summary>

Download both ZIP files using the links above, then install them in a fresh
model/atlas location:

```bash
python - <<'PY'
from platipy.imaging.projects.cardiac.run import install_hybrid_cardiac_from_zip

install_hybrid_cardiac_from_zip(
    "Task400_OPEN_HEART_3d_lowres.zip",
    "open_atlas.zip",
)
PY
```

Keep the installed files available for subsequent runs.

</details>

## Your Images

Use NIfTI volumes (`.nii` or `.nii.gz`). A single-CT segmentation needs only
the image. For dose and longitudinal analysis, create `cohort.json`:

```json
{
  "cases": [
    {
      "case_id": "case001",
      "ct_rt": "data/case001/ct_rt.nii.gz",
      "ct_fu1": "data/case001/ct_fu1.nii.gz",
      "dose": "data/case001/dose.nii.gz",
      "total_dose_gy": 50.0,
      "fractions": 25,
      "dose_type": "effective",
      "modality": "Photon",
      "side": "left"
    }
  ]
}
```

Replace the example paths and treatment values with your own. Add `ct_fu2`
and `ct_fu3` for additional follow-ups. Paths are relative to the registry
file, or can be absolute. See the [full registry example](examples/cohort.example.json).

The dose volume must contain total treatment dose in Gy and be aligned with
the planning CT in physical coordinates. Grid spacing may differ; the workflow
resamples an aligned dose grid but does not register an unaligned volume.

To check the registry and supplied file paths:

```bash
heartdelta registry-validate cohort.json
```

## Running The Workflow

Create `workflow.json` with your registry and output paths:

```json
{
  "registry": "cohort.json",
  "processed_dir": "processed",
  "output_dir": "outputs",
  "segment": true,
  "overwrite": false,
  "alpha_beta": 2.0,
  "hu_normalization": "none"
}
```

To run segmentation, extraction, and the available cohort analyses:

```bash
heartdelta run workflow.json
```

`heartdelta-run workflow.json` is equivalent. Configuration paths are relative
to `workflow.json`. Generated masks are saved under
`processed/<case_id>/<timepoint>/` (`rt`, `fu1`, `fu2`, or `fu3`).

Existing AHA maps are used directly. Missing maps are generated when `segment`
is `true`; cached segmentations are reused unless `overwrite` is `true`.
Result tables and figures are rewritten on each run.

### Reviewed Masks And Wall-Band Measurements

To use reviewed AHA maps, add `aha17_rt`, `aha17_fu1`, and any other timepoint
maps to the registry. Set `"segment": false` to skip automatic segmentation.

Wall-band attenuation measurements also need an LV mask at each timepoint.
After generating and reviewing the masks, create `masks.csv`:

```csv
case_id,timepoint,aha17_path,lv_path
case001,rt,processed/case001/rt/aha17.nii.gz,processed/case001/rt/Ventricle_L.nii.gz
case001,fu1,processed/case001/fu1/aha17.nii.gz,processed/case001/fu1/Ventricle_L.nii.gz
```

Add `"mask_overrides": "masks.csv"` to the workflow configuration and rerun.
Reviewed replacement masks use the same format. Paths are relative to the CSV.
See the [mask override example](examples/mask_overrides.example.csv).

## Parameters And Outputs

The [full workflow example](examples/workflow.example.json) includes additional
options. Adjust its paths, normalization, side filter, and mask overrides
for your cohort before using it.

| Setting | Default | Controls |
| --- | --- | --- |
| `segment` | `true` | Generate missing cardiac masks and AHA maps. |
| `overwrite` | `false` | Regenerate cached masks when segmentation is invoked. |
| `alpha_beta` | `2.0` | α/β in Gy for voxel-wise EQD2. |
| `slab_mm` | `8.0` | Oblique slab thickness for wall-band extraction. |
| `myocardial_band_mm` | `[2.0, 5.0]` | Distance band inside the LV mask. |
| `hu_normalization` | `"none"` | Stored CT intensities or `"trachea_aorta"` calibration for wall-band metrics. |
| `require_complete_three_level_qc` | `true` | Require all three wall-band levels to pass QC. |

EQD2 is calculated voxel-by-voxel before regional averaging. Supply the actual
`fractions` and `total_dose_gy` for each case. For proton cases, set
`"modality": "Proton"` and choose `dose_type` to match the input grid:

| `dose_type` | Handling |
| --- | --- |
| `"physical"` | Multiply the proton grid by `proton_rbe` (default `1.1`). |
| `"effective"` | Use an already RBE-weighted grid without further scaling. |
| `"as_provided"` | Apply no RBE transformation; default when omitted. |

HU calibration requires `trachea_rt`, `aorta_rt`, and corresponding masks for
each follow-up. It maps mean tracheal air to −1000 HU and mean aortic blood to
+50 HU. Calibration details are saved in `native_slab_qc.csv`; segment-wide
metrics use the stored CT intensities.

If the attenuation baseline differs from the planning CT, supply
`attenuation_ct_rt`, `attenuation_aha17_rt`, and `attenuation_lv_rt` for baseline
wall-band measurements. Planning dose still uses `ct_rt` and `aha17_rt`.

Results are saved under `output_dir`:

| Output | Contents |
| --- | --- |
| `segment_metrics.csv`, `level_metrics.csv` | Regional dose, EQD2, gEUD, native HU, and attenuation changes. |
| `aha_qc.csv` | AHA label occupancy and completeness. |
| `processing_status.csv`, `run_summary.json` | Per-case status, errors, and completion counts. |
| `native_slab_metrics.csv`, `native_slab_qc.csv`, `native_slab_level_metrics.csv` | Wall-band measurements, extraction checks, and pooled levels. |
| `attenuation_cohort_qc.csv`, `attenuation_analysis_cohort.csv` | Three-level QC and rows admitted to analysis. |
| `group_tests.csv`, `dose_response.csv`, `mixed_model_terms.csv`, `mixed_model_diagnostics.csv` | Cohort comparisons and dose-response model results. |
| `level_group_boxplots.png`, `level_dose_response.png` | Cohort figures. |

Check `processing_status.csv` and `run_summary.json` after each run. Wall-band
outputs require LV masks; cohort QC requires wall-band data, and model/figure
outputs depend on evaluable observations. Optional tables may be empty.

Supplying `heart_mask`, `lv_mask`, or `lad_mask` also enables whole-structure
DVH measurements in `cardiac_structure_dvh_metrics.csv`. Reference LV/AHA masks
enable segmentation and dose agreement reports.

## Quality Control

AHA labels are 1–6 (basal), 7–12 (mid), 13–16 (apical), and 17 (apex).
The pooled apical level includes 13–17; segment-wide level summaries are
voxel-weighted.

Review cardiac boundaries, orientation, all 17 labels, CT/dose alignment,
dose units, fractionation, and image artifacts before interpreting results.

By default, wall-band cohort analysis requires all three levels to pass:
baseline and follow-up medians must be 0–100 HU, with at least 80% of voxels
in that range. Set `"require_complete_three_level_qc": false` if your protocol
uses a different inclusion rule. Inspect `native_slab_qc.csv` for missing
inputs and failed checks.

## Tests

```bash
pip install -e '.[test]'
python -m pytest tests -q
```

## Citation

Cite the repository URL, version, and commit used. GitHub's **Cite this
repository** control uses [CITATION.cff](CITATION.cff).

## License

[PolyForm Noncommercial License 1.0.0](LICENSE). Commercial use requires a
separate written agreement. External model and atlas licenses apply separately.

For research use. Imaging data and model weights are not bundled.

Maintainer: Patrick Salome
