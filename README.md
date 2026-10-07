# heartdelta-aha17

[![Python](https://img.shields.io/badge/Python-%E2%89%A53.10-3776AB)](pyproject.toml)
[![tests](https://github.com/pgsalome/heartdelta-aha17/actions/workflows/tests.yml/badge.svg)](https://github.com/pgsalome/heartdelta-aha17/actions/workflows/tests.yml)
[![License](https://img.shields.io/badge/use-noncommercial_only-4B5563)](LICENSE)

Cardiac CT segmentation, radiotherapy dose analysis, and longitudinal attenuation
measurements using the American Heart Association 17-segment left-ventricular
model.

Starting from planning and follow-up CT volumes, the workflow can generate
cardiac structure masks and AHA-17 label maps, calculate regional dose and EQD2,
and measure CT attenuation changes. Existing reviewed AHA maps can also be used.

Each CT is analyzed on its native grid with its own AHA map. Planning dose is
sampled on the planning CT; follow-up attenuation is compared by corresponding
AHA segment or level.

## Quick start

### 1. Install

Requires Python 3.10 or newer. Install the segmentation extra to generate
cardiac masks and AHA maps from CT:

```bash
git clone https://github.com/pgsalome/heartdelta-aha17.git
cd heartdelta-aha17
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[segmentation]'
```

Segmentation uses PlatiPy and may download external model weights on first use.
For extraction from existing reviewed maps, use `python -m pip install -e .`.

### 2. Prepare the inputs

Supply CT and dose volumes as NIfTI files (`.nii` or `.nii.gz`). The dose volume
must contain total treatment dose in Gy and share the planning CT's physical
coordinate system. The workflow can resample an aligned dose grid to the CT
spacing; it does not register an unaligned dose volume.

Create `cohort.json` with one entry per case:

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

Replace the paths, fractionation, and case labels with your own data. Add
`ct_fu2` and `ct_fu3` for additional follow-ups. Registry paths are relative to
`cohort.json`, or can be absolute. A fuller example is available in
[`examples/cohort.example.json`](examples/cohort.example.json).

Validate the registry and check that the supplied files exist:

```bash
heartdelta registry-validate cohort.json
```

### 3. Run the workflow

Create `workflow.json`:

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

Then run:

```bash
heartdelta run workflow.json
```

The equivalent standalone command is `heartdelta-run workflow.json`.
Configuration paths are relative to `workflow.json`.

The workflow generates missing AHA maps and extracts segment and level metrics.
Cardiac masks and `aha17.nii.gz` are saved under
`processed/<case_id>/<timepoint>/` (`rt`, `fu1`, `fu2`, or `fu3`).

Supplied AHA maps are used directly. Cached generated masks are reused unless
`"overwrite": true`; this setting applies to scans without a supplied AHA map.
Result tables and figures are rewritten on each run. Check
`outputs/processing_status.csv` and `outputs/run_summary.json` after completion.

## Outputs

Results are written to `output_dir`:

| File | Contents |
| --- | --- |
| `segment_metrics.csv` | Dose, EQD2, gEUD, native HU, and ΔHU for segments 1–17. |
| `level_metrics.csv` | Voxel-weighted basal, mid, and apical summaries. |
| `aha_qc.csv` | Label occupancy and completeness for each AHA map. |
| `processing_status.csv`, `run_summary.json` | Per-case status, errors, and run counts. |
| `native_slab_metrics.csv`, `native_slab_qc.csv` | Wall-band attenuation measurements and extraction/calibration QC. |
| `native_slab_level_metrics.csv` | Pooled wall-band levels joined to planning dose. |
| `attenuation_cohort_qc.csv`, `attenuation_analysis_cohort.csv` | Three-level QC decisions and rows admitted to analysis. |
| `group_tests.csv`, `dose_response.csv` | Proton/photon comparisons and level-specific dose-response fits. |
| `mixed_model_terms.csv`, `mixed_model_diagnostics.csv` | Mixed-effects estimates and fit diagnostics. |
| `level_group_boxplots.png`, `level_dose_response.png` | Cohort attenuation and dose-response figures. |
| `cardiac_structure_dvh_metrics.csv` | DVH metrics when `heart_mask`, `lv_mask`, or `lad_mask` is supplied. |
| `segmentation_validation.csv`, `segmentation_validation_segments.csv`, `segmentation_validation_summary.json` | Agreement results when automatic and reference LV/AHA masks are supplied. |

Wall-band tables require LV masks as described above. QC cohort tables require
wall-band data and `require_complete_three_level_qc`; model and figure outputs
depend on evaluable data. Optional measurements may leave empty CSVs, and some
outputs are written only when their inputs are available.

## Reviewed masks and wall-band extraction

To skip automatic segmentation, add `aha17_rt`, `aha17_fu1`, and any other
available timepoint maps to each registry entry, then set `"segment": false`.
Each map must contain integer labels 1–17 in the corresponding CT's physical
coordinate system.

Myocardial wall-band extraction additionally requires an LV mask at every
analyzed timepoint. To use masks from the first run, create `masks.csv`:

```csv
case_id,timepoint,aha17_path,lv_path
case001,rt,processed/case001/rt/aha17.nii.gz,processed/case001/rt/Ventricle_L.nii.gz
case001,fu1,processed/case001/fu1/aha17.nii.gz,processed/case001/fu1/Ventricle_L.nii.gz
```

Add `"mask_overrides": "masks.csv"` to `workflow.json` and rerun. Review the
masks before using them for analysis. Replacement masks use the same CSV format;
paths resolve from the CSV location. See
[`examples/mask_overrides.example.csv`](examples/mask_overrides.example.csv).

## Analysis options

The full [`workflow example`](examples/workflow.example.json) includes HU
normalization, side filtering, and mask overrides. Adjust its paths and options
for your cohort before using it. Wall-band extraction defaults to an 8-mm slab
and a 2–5 mm band inside the LV mask; change `slab_mm` and
`myocardial_band_mm` to use different settings.

### Dose and fractionation

EQD2 is calculated voxel-by-voxel before regional averaging. Supply the actual
`fractions` and `total_dose_gy` for each case.

For proton cases, set `"modality": "Proton"` and choose `dose_type` to match
the input grid:

| `dose_type` | Dose handling |
| --- | --- |
| `"physical"` | Multiply the proton grid by `proton_rbe` (default `1.1`). |
| `"effective"` | Use an already RBE-weighted grid without further scaling. |
| `"as_provided"` | Apply no RBE transformation. This is the default if omitted. |

### HU normalization

To calibrate wall-band attenuation, set `"hu_normalization": "trachea_aorta"`
and provide `trachea_rt`, `aorta_rt`, and the corresponding masks for each
follow-up (`trachea_fu1`, `aorta_fu1`, etc.). The two-point affine transform maps
mean tracheal air to −1000 HU and mean aortic blood to +50 HU. Calibration details
are recorded in `native_slab_qc.csv`. Segment-wide metrics in
`segment_metrics.csv` use the stored CT intensities.

If the attenuation baseline differs from the planning CT, provide
`attenuation_ct_rt`, `attenuation_aha17_rt`, and `attenuation_lv_rt`. These are
used for baseline wall-band attenuation; `ct_rt` and `aha17_rt` remain the dose
geometry. Baseline calibration masks must match the attenuation CT.

## AHA-17 regions and quality control

| Region | Segments |
| --- | --- |
| Basal | 1–6 |
| Mid | 7–12 |
| Apical | 13–16 |
| Apex | 17 |

The pooled apical level includes segments 13–17. Segment-wide level summaries
are voxel-weighted.

Review cardiac boundaries, basal-to-apical direction, septal/lateral orientation,
all 17 labels, and CT/dose alignment. Check dose units, fractionation, motion,
contrast, reconstruction, and artifacts before interpreting results.

By default, a wall-band follow-up enters cohort analysis only when all three
levels pass: both baseline and follow-up pooled medians must be 0–100 HU, with
at least 80% of voxels in that range. Set
`"require_complete_three_level_qc": false` if your analysis protocol uses a
different inclusion rule. Inspect `native_slab_qc.csv` for missing inputs and
failed extraction checks.

## Tests

```bash
python -m pip install -e '.[test]'
pytest
```

## License and citation

Source-available under the [PolyForm Noncommercial License 1.0.0](LICENSE).
Commercial use requires a separate written agreement.

Cite the repository URL, version, and commit used. Machine-readable citation
metadata are provided in [`CITATION.cff`](CITATION.cff).

For research use. Imaging data and external model weights are not included.

Maintainer: Patrick Salome
