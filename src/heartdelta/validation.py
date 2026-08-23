"""Automatic-versus-reference LV segmentation and dose agreement."""
from __future__ import annotations
import numpy as np
import pandas as pd
import SimpleITK as sitk
from .native import _dose_on_ct, _labels_on_ct
from .metrics import eqd2_values


def dice(a,b):
    a=np.asarray(a)>0; b=np.asarray(b)>0; total=a.sum()+b.sum(); return float(2*np.count_nonzero(a&b)/total) if total else np.nan


def validate_case(case:dict,alpha_beta=2.0) -> tuple[dict,list[dict]]:
    ct_path=case.get("ct_rt") or case.get("baseline_ct"); auto=case.get("lv_mask") or case.get("lv"); reference=case.get("reference_lv_mask") or case.get("manual_lv_mask")
    if not all((ct_path,auto,reference,case.get("dose"))): return {},[]
    ct=sitk.ReadImage(str(ct_path)); a=_labels_on_ct(auto,ct)>0; r=_labels_on_ct(reference,ct)>0; dose=_dose_on_ct(case["dose"],ct); modality=str(case.get("modality",case.get("group",""))).lower()
    if "proton" in modality and str(case.get("dose_type","as_provided")).lower() in {"physical","absorbed"}: dose*=float(case.get("proton_rbe",1.1))
    eqd2=eqd2_values(dose,float(case.get("fractions",case.get("num_fraction_planned",25))),alpha_beta); voxel=float(np.prod(ct.GetSpacing())/1000)
    summary={"case_id":case["case_id"],"lv_dice":dice(a,r),"auto_volume_cc":a.sum()*voxel,"reference_volume_cc":r.sum()*voxel,"auto_mean_dose_eqd2":float(eqd2[a].mean()),"reference_mean_dose_eqd2":float(eqd2[r].mean())}; summary["dose_difference_auto_minus_reference"]=summary["auto_mean_dose_eqd2"]-summary["reference_mean_dose_eqd2"]
    segment_rows=[]; auto_aha=case.get("aha17_rt") or case.get("aha17"); ref_aha=case.get("reference_aha17")
    if auto_aha and ref_aha:
        al=_labels_on_ct(auto_aha,ct); rl=_labels_on_ct(ref_aha,ct)
        for segment in range(1,18): segment_rows.append({"case_id":case["case_id"],"segment":segment,"auto_mean_dose_eqd2":float(eqd2[al==segment].mean()),"reference_mean_dose_eqd2":float(eqd2[rl==segment].mean()),"difference_auto_minus_reference":float(eqd2[al==segment].mean()-eqd2[rl==segment].mean())})
    return summary,segment_rows


def agreement_summary(values:pd.DataFrame) -> dict[str,float]:
    difference=values.dose_difference_auto_minus_reference.dropna(); correlation=values[["auto_mean_dose_eqd2","reference_mean_dose_eqd2"]].corr().iloc[0,1]
    return {"n":len(values),"median_absolute_dose_difference":float(difference.abs().median()),"mean_bias":float(difference.mean()),"lower_95_limit":float(difference.mean()-1.96*difference.std(ddof=1)),"upper_95_limit":float(difference.mean()+1.96*difference.std(ddof=1)),"dose_correlation_r":float(correlation),"dose_correlation_r_squared":float(correlation**2)}
