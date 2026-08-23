import json
from heartdelta.registry import validate


def test_raw_ct_registry_can_be_segmented_later(tmp_path):
    path = tmp_path / "cohort.json"
    path.write_text(json.dumps({"cases":[{"case_id":"x","ct_rt":"ct.nii.gz","dose":"dose.nii.gz"}]}))
    assert validate(path, check_files=False) == []
