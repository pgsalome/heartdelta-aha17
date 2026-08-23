from __future__ import annotations
import numpy as np
import pandas as pd
import SimpleITK as sitk
from .aha import LEVELS,NAMES

def eqd2_factor(total_dose,fractions,alpha_beta=2.0):
    if fractions<=0 or alpha_beta<=0: raise ValueError("fractions and alpha/beta must be positive")
    return (total_dose/fractions+alpha_beta)/(2+alpha_beta)

def eqd2_values(total_voxel_dose,fractions,alpha_beta=2.0):
    """Convert total voxel dose to EQD2 before spatial aggregation."""
    if fractions<=0 or alpha_beta<=0: raise ValueError("fractions and alpha/beta must be positive")
    dose=np.asarray(total_voxel_dose,dtype=float)
    return dose*((dose/fractions+alpha_beta)/(2+alpha_beta))

def geud(values,a):
    x=np.asarray(values,dtype=float); x=x[np.isfinite(x)&(x>0)]
    if not len(x): return np.nan
    if a==0: return float(np.exp(np.mean(np.log(x))))
    return float(np.mean(x**a)**(1/a))

def _read_on_grid(path,reference,nearest=False):
    image=sitk.ReadImage(str(path))
    return sitk.Resample(image,reference,sitk.Transform(),sitk.sitkNearestNeighbor if nearest else sitk.sitkLinear,0,sitk.sitkUInt8 if nearest else sitk.sitkFloat32)

def extract(case,alpha_beta=2.0,a_values=(1,2,3,4,5)):
    baseline=sitk.ReadImage(str(case['baseline_ct'])); segments=_read_on_grid(case['aha17'],baseline,True)
    dose=_read_on_grid(case['dose'],baseline); arrays={'baseline':sitk.GetArrayFromImage(baseline),'dose':sitk.GetArrayFromImage(dose)}
    for tp in (1,2):
        path=case.get(f'followup_{tp}_registered')
        if path: arrays[f'followup_{tp}']=sitk.GetArrayFromImage(_read_on_grid(path,baseline))
    labels=sitk.GetArrayFromImage(segments); total=float(case.get('total_dose_gy',50)); fractions=float(case.get('fractions',total/2)); eqd2=eqd2_values(arrays['dose'],fractions,alpha_beta)
    rows=[]
    for segment in range(1,18):
        mask=labels==segment; n=int(mask.sum()); row={'case_id':case['case_id'],'segment':segment,'segment_name':NAMES[segment],'level':LEVELS[segment],'n_voxels':n}
        for field in ('modality','group','side','target'):
            if field in case: row[field]=case[field]
        if n:
            dv=arrays['dose'][mask]; ev=eqd2[mask]; row.update(mean_dose_gy=float(np.mean(dv)),mean_dose_eqd2_gy=float(np.mean(ev)),baseline_hu_mean=float(np.mean(arrays['baseline'][mask])))
            for a in a_values: row[f'geud_eqd2_a{a:g}']=geud(ev,a)
            for tp in (1,2):
                if f'followup_{tp}' in arrays:
                    value=float(np.mean(arrays[f'followup_{tp}'][mask])); row[f'followup_{tp}_hu_mean']=value; row[f'delta_hu_followup_{tp}']=value-row['baseline_hu_mean']
        rows.append(row)
    return pd.DataFrame(rows)

def aggregate_levels(segments):
    value_cols=[x for x in segments.columns if x.startswith(('mean_dose','baseline_hu','followup_','delta_hu','geud_'))]
    weights=segments.n_voxels.fillna(0).to_numpy(); rows=[]
    for (case,level),group in segments.groupby(['case_id','level'],sort=False):
        idx=group.index; row={'case_id':case,'level':level,'segments':','.join(map(str,group.segment)),'n_voxels':int(group.n_voxels.sum())}
        for field in ('modality','group','side','target'):
            if field in group: row[field]=group[field].iloc[0]
        for column in value_cols:
            valid=group[column].notna()&(group.n_voxels>0); row[column]=float(np.average(group.loc[valid,column],weights=group.loc[valid,'n_voxels'])) if valid.any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)
