# Manuscript method coverage

This inventory maps the draft manuscript to the reusable software. It prevents
an analysis-only convenience path from being mistaken for a complete paper
reproduction.

## Implemented in the clone-to-results workflow

| Manuscript method | Repository implementation |
|---|---|
| Independent cardiac segmentation at planning, FU1, and FU2 | Optional PlatiPy cardiac segmentation in `heartdelta.segmentation`; accepted masks can be supplied instead |
| No CT registration between timepoints | Native-grid extraction in `heartdelta.slab` |
| AHA segments 1–17 and basal/mid/apical-apex pooling | `heartdelta.segmentation`, `heartdelta.aha`, and `heartdelta.native` |
| Voxel-wise cardiac EQD2 with α/β = 2 Gy and patient fraction count | `heartdelta.metrics`; α/β = 2 is the workflow default |
| Full AHA masks for dose | `segment_metrics.csv` and `level_metrics.csv` |
| 8-mm oblique slab at maximum level area | `heartdelta.slab` |
| Band 2–5 mm inward from the outer LV contour | `heartdelta.slab` |
| Trachea/aorta calibration to −1000/+50 HU | `heartdelta.calibration`; coefficients and QC are recorded per timepoint |
| Separate treatment-planning and attenuation baseline CTs | `ct_rt`/`aha17_rt` plus optional `attenuation_*_rt` registry fields |
| Median and ≥80% in-range wall-band QC | `heartdelta.analysis.complete_level_cohort` |
| Same patients at all three levels per follow-up | `attenuation_analysis_cohort.csv` |
| Photon–proton Wilcoxon comparisons | `group_tests.csv` |
| Patient-random-intercept common-slope, dose-by-level, and dose-by-modality models | `mixed_model_terms.csv` and `mixed_model_diagnostics.csv` |
| FU1 photon-only common-slope sensitivity model | `mixed_model_terms.csv` |
| AHA and processing QC audit trails | `aha_qc.csv`, `native_slab_qc.csv`, and `processing_status.csv` |

The oblique-slab implementation was compared with the preserved study outputs:
all 1,455 shared rows matched exactly for voxel count, volume, mean HU, trimmed
mean, median, P95, and fraction of voxels between 0 and 100 HU.

## Inputs or paper-specific products not supplied by this repository

These are not silently reconstructed because doing so would require protected
trial data, clinical judgments, or a study-specific convention:

- DICOM-to-NIfTI and RTSTRUCT conversion. The public workflow starts from
  NIfTI CT, dose, and masks.
- Timepoint-specific trachea and aorta masks. The calibration operation is now
  implemented, but the protected masks are required to run it. The preserved
  canonical slab comparison used stored, untransformed CT voxel values, so a
  normalized rerun is expected to change the archived HU results.
- Physician reference contours and the 53-patient validation dataset. When
  supplied, the software computes Dice, automatic-versus-manual dose agreement,
  limits of agreement, and segment differences; the protected masks are absent.
- Whole-heart, LV, and LAD clinical contours. When supplied, the software
  computes mean dose, V10–V30, and D0.03cc; the protected masks are absent.
- Baseline clinical characteristics and cardiac-adverse-event records used in
  Table 1 and Supplementary Table S1.
- The study's fixed 4.04-Gy photon subgroup threshold. A reusable implementation
  should take this as configuration rather than encode a manuscript result.
- Final journal-layout Figure 1–4 and Supplementary Figure S1–S2 composition.
  The repository generates reusable QC, grouped attenuation, dose-response, and
  bullseye primitives; manuscript panel assembly remains study-specific.
- The manuscript reports R 4.1.2. The public implementation uses Python and
  statsmodels. On the locked 153-row analysis table, it reproduced the draft's
  common, level-specific, modality-interaction, and FU1 photon-only estimates;
  new environments should still verify coefficients and convergence.

## Reproduction contract for this manuscript

A manuscript reproduction registry must identify the 53 left-sided/bilateral
dosimetry cases, their treatment-planning CT/dose/AHA maps, patient-specific
fraction counts, modality, and whole-heart/LV/LAD masks. For attenuation it must
also identify the appropriate free-breathing baseline for photon patients,
native FU1/FU2 CTs, timepoint-specific AHA/LV masks, and the calibrated CT
volumes or calibration provenance. Cohort counts should resolve to 24 photon and
29 proton dosimetry patients, then 14/15 at FU1 and 11/11 at FU2 after locked
three-level QC. A count mismatch is a failed reproduction, not a warning to
ignore.

The organized registry currently yields 26 automated-QC FU1 patients and 21
FU2 patients from its directly referenced masks, whereas the locked manuscript
analysis contains 29 and 22. Before claiming a fresh manuscript reproduction,
the recovered/replaced-mask provenance for those additional three FU1 and one
FU2 cases must be encoded in a versioned registry or mask-override manifest.
The archived 153-row locked analysis table is the source that reproduced the
reported mixed-model coefficients.
