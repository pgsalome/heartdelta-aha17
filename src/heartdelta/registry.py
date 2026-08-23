from __future__ import annotations
import json
from pathlib import Path

PATH_FIELDS=(
    "baseline_ct", "ct_rt", "ct_fu1", "ct_fu2", "ct_fu3", "dose",
    "aha17", "aha17_rt", "aha17_fu1", "aha17_fu2", "aha17_fu3",
    "segment17_platipy", "segment17_platipy_automatic",
    "segment17_platipy_fu1", "segment17_platipy_fu2", "segment17_platipy_fu3",
    "lv_mask", "rv_mask", "la_mask", "heart_mask",
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
