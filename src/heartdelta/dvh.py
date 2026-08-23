"""Whole-structure cardiac DVH metrics."""
from __future__ import annotations
import numpy as np
import pandas as pd
import SimpleITK as sitk
from .metrics import eqd2_values
from .native import _labels_on_ct, _dose_on_ct


def structure_metrics(values:np.ndarray,voxel_cc:float,thresholds=(10,15,20,30)) -> dict[str,float]:
    values=values[np.isfinite(values)]
    if not len(values): return {"mean_dose_eqd2":np.nan,"d0_03cc_eqd2":np.nan,**{f"v{x}_percent":np.nan for x in thresholds}}
    n=max(1,int(np.ceil(.03/voxel_cc)))
    return {"mean_dose_eqd2":float(values.mean()),"d0_03cc_eqd2":float(np.sort(values)[-n:].mean()),**{f"v{x}_percent":float(100*np.mean(values>=x)) for x in thresholds}}


def extract_case(case:dict,alpha_beta=2.0) -> pd.DataFrame:
    ct_path=case.get("ct_rt") or case.get("baseline_ct")
    if not ct_path or not case.get("dose"): return pd.DataFrame()
    ct=sitk.ReadImage(str(ct_path)); dose=_dose_on_ct(case["dose"],ct); modality=str(case.get("modality",case.get("group",""))).lower(); dose_type=str(case.get("dose_type","as_provided")).lower()
    if "proton" in modality and dose_type in {"physical","absorbed"}: dose*=float(case.get("proton_rbe",1.1))
    fractions=float(case.get("fractions",case.get("num_fraction_planned",25))); eqd2=eqd2_values(dose,fractions,alpha_beta); voxel_cc=float(np.prod(ct.GetSpacing())/1000); rows=[]
    candidates={"heart":case.get("heart_mask") or case.get("heart"),"lv":case.get("lv_mask") or case.get("lv"),"lad":case.get("lad_mask") or case.get("lad")}
    for structure,path in candidates.items():
        if not path: continue
        mask=_labels_on_ct(path,ct)>0; row={"case_id":case["case_id"],"structure":structure,"volume_cc":float(mask.sum()*voxel_cc),**structure_metrics(eqd2[mask],voxel_cc)}
        for field in ("modality","group","side"): row[field]=case.get(field,"")
        rows.append(row)
    return pd.DataFrame(rows)
