"""Two-point CT HU calibration using tracheal air and aortic blood."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import SimpleITK as sitk
from .native import _labels_on_ct


def parameters(ct:sitk.Image,trachea_path,aorta_path,target_air=-1000.0,target_blood=50.0,min_voxels=50):
    values=sitk.GetArrayFromImage(sitk.Cast(ct,sitk.sitkFloat32)).astype(float)
    trachea=_labels_on_ct(trachea_path,ct)>0; aorta=_labels_on_ct(aorta_path,ct)>0
    air=values[trachea&np.isfinite(values)]; blood=values[aorta&np.isfinite(values)]
    if len(air)<min_voxels or len(blood)<min_voxels: raise ValueError(f"insufficient calibration voxels: trachea={len(air)}, aorta={len(blood)}")
    mean_air=float(air.mean()); mean_blood=float(blood.mean()); denominator=mean_air-mean_blood
    if not np.isfinite(denominator) or abs(denominator)<100: raise ValueError(f"invalid trachea/aorta separation: {denominator:.3f} HU")
    scale=(target_air-target_blood)/denominator; intercept=target_blood-mean_blood*scale
    return {"mean_trachea_hu":mean_air,"mean_aorta_hu":mean_blood,"target_air_hu":float(target_air),"target_blood_hu":float(target_blood),"scale":float(scale),"intercept":float(intercept),"trachea_voxels":int(len(air)),"aorta_voxels":int(len(blood))}


def apply(values,calibration): return np.asarray(values,dtype=float)*calibration["scale"]+calibration["intercept"]


def normalize_image(ct_path,trachea_path,aorta_path,output_path,target_air=-1000.0,target_blood=50.0):
    ct=sitk.ReadImage(str(ct_path)); calibration=parameters(ct,trachea_path,aorta_path,target_air,target_blood); array=apply(sitk.GetArrayFromImage(ct),calibration).astype(np.float32); output=sitk.GetImageFromArray(array); output.CopyInformation(ct); target=Path(output_path); target.parent.mkdir(parents=True,exist_ok=True); sitk.WriteImage(output,str(target)); target.with_suffix(target.suffix+".calibration.json").write_text(json.dumps(calibration,indent=2)+"\n"); return calibration


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--ct",required=True); parser.add_argument("--trachea",required=True); parser.add_argument("--aorta",required=True); parser.add_argument("--output",required=True); parser.add_argument("--target-air",type=float,default=-1000); parser.add_argument("--target-blood",type=float,default=50); args=parser.parse_args(argv); print(json.dumps(normalize_image(args.ct,args.trachea,args.aorta,args.output,args.target_air,args.target_blood),indent=2)); return 0
