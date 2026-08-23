from __future__ import annotations
import json
from pathlib import Path

PATH_FIELDS=("baseline_ct","dose","aha17","lv_mask","rv_mask","followup_1_registered","followup_2_registered")

def load(path):
    source=Path(path); payload=json.loads(source.read_text()); cases=payload.get("cases") if isinstance(payload,dict) else payload
    if not isinstance(cases,list): raise ValueError("Registry must be a list or contain a 'cases' list")
    out=[]
    for item in cases:
        case=dict(item)
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
        if not case.get("aha17") and not (case.get("lv_mask") and case.get("rv_mask")):
            errors.append(f"{cid}: provide aha17 or both lv_mask and rv_mask")
        if check_files:
            for field in PATH_FIELDS:
                value=case.get(field)
                if value and not Path(value).is_file(): errors.append(f"{cid}: {field} does not exist: {value}")
    return errors

