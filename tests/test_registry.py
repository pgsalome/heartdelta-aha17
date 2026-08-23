import json
from heartdelta.registry import apply_mask_overrides, load, validate


def test_raw_ct_registry_can_be_segmented_later(tmp_path):
    path = tmp_path / "cohort.json"
    path.write_text(json.dumps({"cases":[{"case_id":"x","ct_rt":"ct.nii.gz","dose":"dose.nii.gz"}]}))
    assert validate(path, check_files=False) == []


def test_mask_overrides_replace_timepoint_paths(tmp_path):
    registry = tmp_path / "cohort.json"
    registry.write_text(json.dumps({"cases":[{"case_id":"x","ct_rt":"ct.nii.gz","dose":"dose.nii.gz"}]}))
    overrides = tmp_path / "overrides.csv"
    overrides.write_text(
        "case_id,timepoint,aha17_path,lv_path,source_method\n"
        "x,fu1,masks/aha.nii.gz,masks/lv.nii.gz,reviewed\n"
    )
    cases = apply_mask_overrides(load(registry), overrides)
    assert cases[0]["aha17_fu1"] == str((tmp_path / "masks/aha.nii.gz").resolve())
    assert cases[0]["lv_fu1_auto"] == str((tmp_path / "masks/lv.nii.gz").resolve())
    assert cases[0]["mask_override_provenance"]["fu1"]["source_method"] == "reviewed"
