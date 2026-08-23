# heartdelta-aha17

[![Python](https://img.shields.io/badge/Python-%E2%89%A53.10-3776AB)](pyproject.toml)
[![License](https://img.shields.io/badge/use-noncommercial_only-4B5563)](LICENSE)

## Put longitudinal cardiac CT and dose on the same 17-segment map

Whole-heart and whole-LV averages can hide spatially concentrated treatment
effects. `heartdelta-aha17` organizes planning dose and longitudinal CT
attenuation using the American Heart Association 17-segment left-ventricular
model.

The package handles two related jobs:

1. Generate an AHA-17 label volume from LV and RV masks using an anatomical
   long axis and an RV-facing angular landmark.
2. Extract segment-level and basal/mid/apical summaries from an existing or
   newly generated AHA map.

Automatic cardiac segmentation is deliberately not tied to one model. Bring
physician contours, TotalSegmentator/Platipy outputs, or another reviewed LV/RV
pair; this toolkit owns the mapping and measurement steps that follow.

## Installation

```bash
git clone https://github.com/pgsalome/heartdelta-aha17.git
cd heartdelta-aha17
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Cohort registry

One JSON record represents one patient or independent analysis case:

```json
{
  "cases": [
    {
      "case_id": "case001",
      "baseline_ct": "/cohort/case001/planning_ct.nii.gz",
      "followup_1_registered": "/cohort/case001/ct_6m_on_planning.nii.gz",
      "followup_2_registered": "/cohort/case001/ct_12m_on_planning.nii.gz",
      "dose": "/cohort/case001/dose_on_planning.nii.gz",
      "aha17": "/cohort/case001/aha17.nii.gz",
      "fractions": 25,
      "total_dose_gy": 50,
      "modality": "Proton",
      "side": "left"
    }
  ]
}
```

Instead of `aha17`, provide both `lv_mask` and `rv_mask` when the segment map
needs to be generated. Paths may be absolute or relative to the registry file.

```bash
heartdelta registry-validate my_heart_cohort.json
```

## Generate the AHA geometry

```bash
heartdelta aha-generate \
  --registry my_heart_cohort.json \
  --shell-mm 8 \
  --output results/aha_maps
```

The generator:

- constructs an inner LV shell in physical millimetres;
- estimates the LV long axis with PCA;
- uses the RV-facing surface as an angular reference;
- divides basal and mid levels into six sectors each;
- divides the apical level into four sectors;
- assigns the terminal apex to segment 17.

This is a geometric approximation, not a substitute for contour review.
Inspect basal/apical direction, septal/lateral orientation, segment occupancy,
and agreement with cardiac anatomy before accepting a generated map.

Update the registry to point `aha17` at each accepted output before extraction.

## Extract longitudinal measurements

```bash
heartdelta extract \
  --registry my_heart_cohort.json \
  --alpha-beta 3 \
  --output results/segment_analysis
```

The command resamples dose and registered follow-ups onto the planning CT grid
and writes:

- `segment_metrics.csv`: one row per case and AHA segment;
- `level_metrics.csv`: voxel-weighted basal, mid, and apical/apex summaries.

Each segment row contains voxel count, mean planned dose, mean EQD2, baseline
HU, follow-up HU, ΔHU, and gEUD values for several `a` parameters. Cohort labels
from the registry are carried into both tables for downstream models.

EQD2 uses the registry's total dose and fraction count. If those fields are
absent, the extractor assumes 50 Gy in 25 fractions; provide the real values for
any definitive analysis.

## Segment numbering

```text
Basal:   1 anterior, 2 anteroseptal, 3 inferoseptal,
         4 inferior, 5 inferolateral, 6 anterolateral

Mid:     7 anterior, 8 anteroseptal, 9 inferoseptal,
        10 inferior, 11 inferolateral, 12 anterolateral

Apical: 13 anterior, 14 septal, 15 inferior, 16 lateral
Apex:   17
```

The level table combines segments 1–6, 7–12, and 13–17. Its averages are
weighted by the number of contributing voxels, rather than giving a tiny apex
segment the same influence as an entire basal sector.

## What this repository does not assume

- It does not require a particular cardiac autosegmentation network.
- It does not perform raw DICOM import.
- It does not silently register follow-up CT; inputs named `*_registered` must
  already be aligned and QC-approved.
- It does not prescribe a statistical model for your cohort.
- It does not include patient data, manuscript figures, or unpublished results.

Keeping these boundaries explicit makes it easier to substitute segmentation
or registration methods without changing the downstream AHA measurements.

## Quality checks that matter

- Confirm that CT, dose, AHA labels, and follow-ups describe the same physical
  coordinate system.
- Require all labels 1–17 for a complete map and inspect unexpected empty or
  very small segments.
- View the RV-facing reference and septal/lateral assignment on several slices.
- Inspect longitudinal registration specifically along the myocardium.
- Report whether measurements use a synthetic LV shell, a myocardium contour,
  or another structure definition.
- Keep patient-level independence in subsequent regression or mixed models.

## Tests

```bash
pip install -e '.[test]'
pytest
```

The tests cover all 17 label assignments, EQD2 behavior, and gEUD calculation.
Real cardiac geometry still requires visual and clinical QC.

## Research and privacy notice

This software is for research use. It is not a clinical segmentation system,
treatment-planning product, or diagnostic device. No patient data should be
committed to the repository; `data/`, `results/`, and `outputs/` are ignored.

There is not yet an associated paper for this standalone toolkit. Cite the
repository release and commit used until a formal software citation is issued.

## License

This software is source-available under the
[PolyForm Noncommercial License 1.0.0](LICENSE).

Noncommercial use, modification, and redistribution are permitted subject
to the license terms. Commercial use requires a separate written agreement.

This is not an OSI-approved open-source license.

## Maintainer

Patrick Salome
