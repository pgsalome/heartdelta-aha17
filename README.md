# heartdelta-aha17

[![Python](https://img.shields.io/badge/Python-%E2%89%A53.10-3776AB)](pyproject.toml)
[![tests](https://github.com/pgsalome/heartdelta-aha17/actions/workflows/tests.yml/badge.svg)](https://github.com/pgsalome/heartdelta-aha17/actions/workflows/tests.yml)
[![License](https://img.shields.io/badge/use-noncommercial_only-4B5563)](LICENSE)

Reproducible longitudinal cardiac CT and radiotherapy-dose analysis using the
American Heart Association 17-segment left-ventricular model.

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
- `group_tests.csv` and `dose_response.csv`;
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

## Maintainer

Patrick Salome
