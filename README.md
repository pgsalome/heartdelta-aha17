# heartdelta-aha17

[![Python](https://img.shields.io/badge/Python-%E2%89%A53.10-3776AB)](pyproject.toml)
[![tests](https://github.com/pgsalome/heartdelta-aha17/actions/workflows/tests.yml/badge.svg)](https://github.com/pgsalome/heartdelta-aha17/actions/workflows/tests.yml)
[![License](https://img.shields.io/badge/use-noncommercial_only-4B5563)](LICENSE)

Reproducible longitudinal cardiac CT and radiotherapy-dose analysis using the
American Heart Association 17-segment left-ventricular model.

The workflow accompanies the manuscript *Left Ventricular Dose Distribution
After Proton Versus Photon Regional Nodal Irradiation for Breast Cancer: AHA
17-Segment Dosimetry*. The code is cohort-reusable; study data are not included.
See [MANUSCRIPT_COVERAGE.md](MANUSCRIPT_COVERAGE.md) for an exact method-to-code
inventory and the remaining manuscript-specific inputs.

## Manuscript reproduction status

The core computations are reproducible: voxel-wise cardiac EQD2, 17-segment
dosimetry, native 8-mm/2–5-mm wall-band extraction, complete three-level QC,
Wilcoxon comparisons, and the manuscript mixed-effects models. The locked
153-row attenuation/model table reproduces the reported coefficients.

A completely fresh paper reproduction is not yet claimable from the organized
registry alone. Two provenance items must be resolved first:

1. Trachea/aorta normalization is now implemented, but the preserved canonical
   slab values match the untransformed stored CT voxels exactly. A fresh
   normalized rerun will therefore be a method correction, not a byte-for-byte
   reproduction of those archived values.
2. Directly referenced masks yield 26 FU1 and 21 FU2 complete-QC patients,
   whereas the locked manuscript cohorts contain 29 and 22. Recovered/replaced
   masks for those additional cases need a versioned override manifest.

The repository does not hide these count mismatches or substitute archived CSVs
for raw-image processing.

This repository is intended for a new cohort, not only for regenerating one
figure. Starting from planning and follow-up CT volumes, it can create cardiac
chamber masks and native-grid AHA-17 maps, extract planning dose and longitudinal
attenuation, run cohort comparisons, and generate level-based figures.

## Analysis design

Follow-up CT intensities remain on their native grids. Each timepoint receives
its own reviewed cardiac/AHA segmentation, and corresponding AHA regions are
compared longitudinally. Planning dose is sampled on the planning CT AHA map.
This avoids interpolating follow-up HU through a deformable warp.

```text
planning CT ── cardiac segmentation ── AHA-17 ── planning dose metrics
follow-up CT ─ cardiac segmentation ── AHA-17 ── native HU metrics
                                      │
                                      └─ segment/level changes ─ statistics + figures
```

## Installation

For extraction from existing reviewed AHA maps:

```bash
git clone https://github.com/pgsalome/heartdelta-aha17.git
cd heartdelta-aha17
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

To generate cardiac structures and AHA maps automatically, install the PlatiPy
cardiac dependencies:

```bash
pip install -e '.[segmentation]'
```

Segmentation models may download external weights on first use. A GPU is useful
but is not required for downstream extraction.

## Input registry

Copy `examples/cohort.example.json`. Each case may contain:

- `ct_rt`: planning CT NIfTI;
- `ct_fu1`, `ct_fu2`, `ct_fu3`: available follow-up CTs;
- `dose`: physical dose on, or relatable to, the planning grid;
- `total_dose_gy` and `fractions`: used to calculate EQD2;
- labels such as `modality`, `side`, `breathing`, and `delivery_technique`.

If accepted AHA maps exist, provide `aha17_rt`, `aha17_fu1`, etc. The pipeline
then skips segmentation. Corrected masks can therefore replace automatic masks
without changing downstream code. Relative paths resolve from the registry.

Cardiac EQD2 defaults to α/β = 2 Gy and is calculated voxel-by-voxel before
regional averaging. Set `dose_type` to `physical` for a proton grid that still
requires `proton_rbe` (1.1 in the manuscript), or to `effective` when the grid
is already RBE-weighted. `as_provided` applies no RBE transformation.

When the treatment-planning CT and attenuation baseline differ—for example,
breath-hold dosimetry but free-breathing longitudinal comparison—keep `ct_rt`
and `aha17_rt` as the dose geometry and additionally provide
`attenuation_ct_rt`, `attenuation_aha17_rt`, and `attenuation_lv_rt`.

## Complete workflow

Copy and edit `examples/workflow.example.json`, then run:

```bash
heartdelta-run my_workflow.json
# equivalent
heartdelta run my_workflow.json
```

The workflow is restartable. Products are reused unless `"overwrite": true`.
Set `"segment": false` for an analysis-only run with reviewed AHA maps.

It writes:

- `segment_metrics.csv`: dose, EQD2, gEUD, native HU, and ΔHU for segments 1–17;
- `level_metrics.csv`: voxel-weighted basal, mid, and apical/apex summaries;
- `native_slab_metrics.csv`: validated 8-mm oblique-slab measurements within
  the configurable 2–5 mm myocardial band;
- `native_slab_level_metrics.csv`: pooled basal, mid, and apical values joined
  to planning dose and used by the supplied cohort statistics and figures;
- `aha_qc.csv`: label occupancy and completeness;
- `attenuation_cohort_qc.csv` and `attenuation_analysis_cohort.csv`: complete
  three-level QC and the exact rows admitted to modeling;
- `group_tests.csv`, `dose_response.csv`, `mixed_model_terms.csv`, and
  `mixed_model_diagnostics.csv`;
- `cardiac_structure_dvh_metrics.csv` when heart/LV/LAD masks are supplied;
- `segmentation_validation.csv`, segment differences, and an agreement summary
  when physician/reference LV and AHA masks are supplied;
- `level_group_boxplots.png` and `level_dose_response.png`;
- `run_summary.json` with completion counts.

Validate configuration inputs separately with:

```bash
heartdelta registry-validate cohort.json
```

The original lightweight registered-image interface remains available:

```bash
heartdelta aha-generate --registry cohort.json --output processed/aha_maps
heartdelta extract --registry cohort.json --output outputs/registered_analysis
```

For new longitudinal cohorts, use the native-grid workflow.

## AHA numbering

```text
Basal:   1 anterior, 2 anteroseptal, 3 inferoseptal,
         4 inferior, 5 inferolateral, 6 anterolateral
Mid:     7 anterior, 8 anteroseptal, 9 inferoseptal,
        10 inferior, 11 inferolateral, 12 anterolateral
Apical: 13 anterior, 14 septal, 15 inferior, 16 lateral
Apex:   17
```

Level summaries are voxel-weighted, so a small apex does not receive the same
influence as an entire basal sector.

## Quality control

Inspect cardiac boundaries, basal-to-apical direction, septal/lateral
orientation, occupancy of all 17 labels, CT/dose physical-coordinate agreement,
dose units, fractionation, motion, contrast, reconstruction, and metal artifact.
The automated QC table does not replace visual clinical review.

The manuscript excludes a follow-up unless both baseline and follow-up pooled
wall-band medians are 0–100 HU, at least 80% of their voxels are within 0–100
HU, and all three levels pass. This is the default workflow behavior. For a
different prospective protocol it can be disabled explicitly with
`"require_complete_three_level_qc": false`.

HU calibration is deliberately not guessed. If scanner-reference calibration
is required, supply timepoint-specific `trachea_*` and `aorta_*` masks and set
`"hu_normalization": "trachea_aorta"`. The pipeline uses the historical
two-point affine transform based on full-mask means, mapping tracheal air to
−1000 HU and aortic blood to +50 HU. Missing masks, too few voxels, or inadequate
landmark separation fail QC; the measured means, scale, intercept, and voxel
counts are written to `native_slab_qc.csv`.

To create a normalized CT independently:

```bash
heartdelta-normalize \
  --ct ct.nii.gz --trachea trachea.nii.gz --aorta aorta.nii.gz \
  --output ct_normalized.nii.gz
```

This writes the normalized volume plus a calibration JSON sidecar.

This software is for research use and is not a treatment-planning or diagnostic
device. Patient data and local paths must not be committed.

## Testing

```bash
pip install -e '.[test]'
pytest
```

Tests cover AHA assignments, EQD2/gEUD, registry behavior, and voxel-weighted
aggregation. Run an end-to-end acceptance case locally because imaging data and
model weights are not distributed here.

## License

This software is source-available under the
[PolyForm Noncommercial License 1.0.0](LICENSE).

Noncommercial use, modification, and redistribution are permitted subject to
the license terms. Commercial use requires a separate written agreement.

This is not an OSI-approved open-source license.

## Citation

Until the associated manuscript has a final citation, cite the repository URL,
release/version, and commit used. Machine-readable metadata are provided in
[`CITATION.cff`](CITATION.cff).

## Maintainer

Patrick Salome
