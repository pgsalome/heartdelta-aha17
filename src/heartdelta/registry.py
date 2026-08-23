from __future__ import annotations
import json
import csv
import os
from pathlib import Path

PATH_FIELDS=(
    "baseline_ct", "ct_rt", "ct_fu1", "ct_fu2", "ct_fu3", "dose",
    "attenuation_ct_rt", "attenuation_aha17_rt", "attenuation_lv_rt",
    "aha17", "aha17_rt", "aha17_fu1", "aha17_fu2", "aha17_fu3",
    "segment17_platipy", "segment17_platipy_automatic",
    "segment17_platipy_fu1", "segment17_platipy_fu2", "segment17_platipy_fu3",
    "lv_mask", "rv_mask", "la_mask", "heart_mask",
    "heart", "lv", "lad", "lad_mask", "reference_lv_mask", "manual_lv_mask", "reference_aha17",
    "trachea_rt", "aorta_rt", "trachea_fu1", "aorta_fu1", "trachea_fu2", "aorta_fu2", "trachea_fu3", "aorta_fu3", "trachea_auto", "aorta_auto",
    "followup_1_registered", "followup_2_registered",
)

def load(path):
    source=Path(path); payload=json.loads(source.read_text()); cases=payload.get("cases") if isinstance(payload,dict) else payload
    if not isinstance(cases,list): raise ValueError("Registry must be a list or contain a 'cases' list")
    out=[]
    for item in cases:
        case=dict(item)
        case.setdefault("case_id", case.get("patient_id"))
        case.setdefault("baseline_ct", case.get("ct_rt"))
        case.setdefault("aha17_rt", case.get("segment17_platipy_automatic") or case.get("segment17_platipy"))
        for tp in (1, 2, 3):
            case.setdefault(f"aha17_fu{tp}", case.get(f"segment17_platipy_fu{tp}"))
        for key in PATH_FIELDS:
            value=case.get(key)
            if value and not Path(value).is_absolute(): case[key]=str((source.parent/value).resolve())
        out.append(case)
    return out


def apply_mask_overrides(cases, path):
    """Apply an audited CSV of accepted timepoint-specific AHA/LV masks.

    Required columns are ``case_id`` (or ``patient_id``), ``timepoint`` (or
    ``phase``), and ``aha17_path``. ``lv_path`` is optional. Relative paths
    resolve from the CSV and environment variables are expanded.
    """
    source = Path(path)
    by_id = {str(case.get("case_id") or case.get("patient_id")): case for case in cases}
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row_number, row in enumerate(rows, start=2):
        case_id = str(row.get("case_id") or row.get("patient_id") or "").strip()
        timepoint = str(row.get("timepoint") or row.get("phase") or "").strip().lower()
        if case_id not in by_id:
            raise ValueError(f"{source}:{row_number}: unknown case_id {case_id!r}")
        if timepoint not in {"rt", "fu1", "fu2", "fu3"}:
            raise ValueError(f"{source}:{row_number}: invalid timepoint {timepoint!r}")
        aha17 = str(row.get("aha17_path") or "").strip()
        if not aha17:
            raise ValueError(f"{source}:{row_number}: missing aha17_path")

        def resolve(value):
            expanded = Path(os.path.expandvars(os.path.expanduser(value)))
            return str(expanded if expanded.is_absolute() else (source.parent / expanded).resolve())

        case = by_id[case_id]
        case[f"aha17_{timepoint}"] = resolve(aha17)
        lv = str(row.get("lv_path") or "").strip()
        if lv:
            case["attenuation_lv_rt" if timepoint == "rt" else f"lv_{timepoint}_auto"] = resolve(lv)
        case.setdefault("mask_override_provenance", {})[timepoint] = {
            key: value for key, value in row.items() if value not in (None, "")
        }
    return cases

def validate(path,check_files=True):
    errors=[]; seen=set()
    for case in load(path):
        cid=str(case.get("case_id",""))
        if not cid: errors.append("case without case_id")
        if cid in seen: errors.append(f"{cid}: duplicate case_id")
        seen.add(cid)
        for field in ("baseline_ct","dose"):
            if not case.get(field): errors.append(f"{cid}: missing {field}")
        # An AHA map is optional at validation time: the complete workflow can
        # generate cardiac structures and AHA labels directly from the CT.
        if check_files:
            for field in PATH_FIELDS:
                value=case.get(field)
                if value and not Path(value).is_file(): errors.append(f"{cid}: {field} does not exist: {value}")
    return errors
