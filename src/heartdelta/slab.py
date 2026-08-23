"""Oblique native-CT myocardial slab extraction used by the validated workflow."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import SimpleITK as sitk
from scipy.ndimage import distance_transform_edt
from .native import _labels_on_ct
from .calibration import apply as apply_calibration, parameters as calibration_parameters

LEVEL_SEGMENTS={"basal":tuple(range(1,7)),"mid":tuple(range(7,13)),"apical":tuple(range(13,18))}


def _physical(indices_zyx, image, start):
    xyz=(indices_zyx.astype(float)+start)[:,::-1]; scaled=xyz*np.asarray(image.GetSpacing()); direction=np.asarray(image.GetDirection()).reshape(3,3)
    return np.asarray(image.GetOrigin())+scaled@direction.T


def _centroid(labels, selected, image, start):
    indices=np.argwhere(np.isin(labels,selected)); return _physical(indices,image,start).mean(0) if len(indices) else None


def _axis(labels,image,start):
    basal=_centroid(labels,range(1,7),image,start); apical=_centroid(labels,range(13,17),image,start)
    if apical is None: apical=_centroid(labels,(17,),image,start)
    if basal is not None and apical is not None and np.linalg.norm(basal-apical)>1: return (basal-apical)/np.linalg.norm(basal-apical),"anatomical"
    points=_physical(np.argwhere(labels>0),image,start); _,_,vectors=np.linalg.svd(points-points.mean(0),full_matrices=False); return vectors[-1]/np.linalg.norm(vectors[-1]),"pca_fallback"


def _center(projections,bin_mm=1.0):
    low=np.floor(projections.min()/bin_mm)*bin_mm; high=np.ceil(projections.max()/bin_mm)*bin_mm+bin_mm; counts,edges=np.histogram(projections,bins=np.arange(low,high+bin_mm*.5,bin_mm)); smooth=np.convolve(counts.astype(float),np.ones(3),mode="same"); index=int(np.median(np.flatnonzero(smooth==smooth.max()))); return float((edges[index]+edges[index+1])/2)


def _trimmed(values,fraction=.1):
    values=np.sort(values); trim=int(np.floor(len(values)*fraction)); return float(values[trim:-trim].mean()) if trim and 2*trim<len(values) else float(values.mean())


def extract_timepoint(case:dict,timepoint:str,level:str,slab_mm:float=8.0,band=(2.0,5.0),normalization="none",target_air=-1000.0,target_blood=50.0):
    ct_path=case.get(f"ct_{timepoint}") if timepoint!="rt" else case.get("attenuation_ct_rt") or case.get("ct_rt") or case.get("baseline_ct")
    aha_path=case.get(f"aha17_{timepoint}") if timepoint!="rt" else case.get("attenuation_aha17_rt") or case.get("aha17_rt") or case.get("aha17")
    lv_path=case.get(f"lv_{timepoint}_auto") if timepoint!="rt" else case.get("attenuation_lv_rt") or case.get("lv_auto") or case.get("lv") or case.get("lv_mask")
    if not all(x and Path(x).is_file() for x in (ct_path,aha_path,lv_path)): return [],{"status":"missing_input"}
    ct_image=sitk.ReadImage(str(ct_path)); ct=sitk.GetArrayFromImage(sitk.Cast(ct_image,sitk.sitkFloat32)); calibration={}
    if normalization=="trachea_aorta":
        trachea=case.get(f"trachea_{timepoint}") or (case.get("trachea_auto") if timepoint=="rt" else None); aorta=case.get(f"aorta_{timepoint}") or (case.get("aorta_auto") if timepoint=="rt" else None)
        if not trachea or not aorta: raise ValueError(f"missing {timepoint} trachea/aorta calibration masks")
        calibration=calibration_parameters(ct_image,trachea,aorta,target_air,target_blood); ct=apply_calibration(ct,calibration)
    labels=np.rint(_labels_on_ct(aha_path,ct_image)).astype(np.int16); lv=_labels_on_ct(lv_path,ct_image)>0
    selected=LEVEL_SEGMENTS[level]
    if not set(selected).issubset(set(np.unique(labels))): raise ValueError(f"incomplete {level} labels")
    overlap=np.count_nonzero((labels>0)&lv)/max(np.count_nonzero(labels>0),1)
    if overlap<.8: raise ValueError(f"AHA/LV overlap {overlap:.3f}")
    loc=np.where(lv|(labels>0)); crop=tuple(slice(max(0,int(x.min())-2),min(lv.shape[i],int(x.max())+3)) for i,x in enumerate(loc)); start=np.asarray([x.start for x in crop],float); cct=ct[crop]; clabels=labels[crop]; clv=lv[crop]
    axis,method=_axis(clabels,ct_image,start); ring=np.argwhere(np.isin(clabels,selected)); center=_center(_physical(ring,ct_image,start)@axis); allidx=np.indices(clabels.shape).reshape(3,-1).T; slab=(np.abs(_physical(allidx,ct_image,start)@axis-center)<=slab_mm/2).reshape(clabels.shape)
    distance=distance_transform_edt(clv,sampling=np.asarray(ct_image.GetSpacing()[::-1])); common=slab&clv&(distance>=band[0])&(distance<=band[1]); voxel=float(np.prod(ct_image.GetSpacing())); rows=[]
    for segment in (*selected,0):
        mask=common&(np.isin(clabels,selected) if segment==0 else clabels==segment); values=cct[mask]; values=values[np.isfinite(values)]
        rows.append({"case_id":str(case["case_id"]),"timepoint":timepoint,"level":level,"segment":segment if segment else pd.NA,"region":f"segments{selected[0]}-{selected[-1]}_pooled" if not segment else f"segment{segment}","n_voxels":len(values),"volume_cc":len(values)*voxel/1000,"hu_mean":float(values.mean()) if len(values) else np.nan,"hu_trimmed_mean":_trimmed(values) if len(values) else np.nan,"hu_median":float(np.median(values)) if len(values) else np.nan,"hu_p95":float(np.percentile(values,95)) if len(values) else np.nan,"fraction_hu_0_100":float(np.mean((values>=0)&(values<=100))) if len(values) else np.nan})
    return rows,{"status":"ok","hu_normalization":normalization,"aha_lv_overlap_fraction":overlap,"axis_method":method,"slab_center_projection_mm":center,"ring_voxels":len(ring),"slab_mm":slab_mm,"band_low_mm":band[0],"band_high_mm":band[1],**calibration}


def extract_case(case:dict,slab_mm=8.0,band=(2.0,5.0),normalization="none",target_air=-1000.0,target_blood=50.0):
    rows=[]; qc=[]
    for tp in ("rt","fu1","fu2","fu3"):
        for level in LEVEL_SEGMENTS:
            try: part,status=extract_timepoint(case,tp,level,slab_mm,band,normalization,target_air,target_blood)
            except (RuntimeError,ValueError) as exc: part,status=[],{"status":"failed_qc","error":str(exc)}
            rows.extend(part); qc.append({"case_id":case["case_id"],"timepoint":tp,"level":level,**status})
    frame=pd.DataFrame(rows)
    if not frame.empty:
        keys=["level","region","segment"]; base_values=["hu_mean","hu_trimmed_mean","hu_median","hu_p95","fraction_hu_0_100"]; base=frame[frame.timepoint=="rt"][keys+base_values].rename(columns={x:f"baseline_{x}" for x in base_values}); frame=frame.merge(base,on=keys,how="left")
        for x in ("hu_mean","hu_trimmed_mean","hu_median","hu_p95"): frame[f"delta_{x}"]=frame[x]-frame[f"baseline_{x}"]
        for field in ("group","modality","side","breathing","delivery_technique","target"):
            if field in case: frame[field]=case[field]
    return frame,pd.DataFrame(qc)
