from __future__ import annotations
import math
import numpy as np
import SimpleITK as sitk
from scipy.spatial import cKDTree

NAMES={1:"Basal anterior",2:"Basal anteroseptal",3:"Basal inferoseptal",4:"Basal inferior",5:"Basal inferolateral",6:"Basal anterolateral",7:"Mid anterior",8:"Mid anteroseptal",9:"Mid inferoseptal",10:"Mid inferior",11:"Mid inferolateral",12:"Mid anterolateral",13:"Apical anterior",14:"Apical septal",15:"Apical inferior",16:"Apical lateral",17:"Apex"}
LEVELS={**{x:"basal" for x in range(1,7)},**{x:"mid" for x in range(7,13)},**{x:"apical" for x in range(13,18)}}

def assign_segments(z,angle):
    z=np.asarray(z); angle=np.asarray(angle)%360; low=z.min(); span=z.max()-low
    if span<=0: raise ValueError("LV long-axis extent is zero")
    labels=np.zeros(len(z),np.uint8); basal=low+2*span/3; mid=low+span/3; apex=low+.1*span
    six=np.floor(((angle+30)%360)/60).astype(int)
    labels[z>=basal]=np.take([1,2,3,4,5,6],six[z>=basal])
    use=(z>=mid)&(z<basal); labels[use]=np.take([7,8,9,10,11,12],six[use])
    use=(z>=apex)&(z<mid); labels[use]=13+np.floor(((angle[use]+15)%360)/90).astype(int)
    labels[z<apex]=17
    return labels

def _points(image,array):
    idx=np.argwhere(array); xyz=idx[:,[2,1,0]].astype(float)
    spacing=np.asarray(image.GetSpacing()); origin=np.asarray(image.GetOrigin()); direction=np.asarray(image.GetDirection()).reshape(3,3)
    return idx,origin+(xyz*spacing)@direction.T

def generate(lv_path,rv_path,output,shell_thickness_mm=8.0):
    lv=sitk.ReadImage(str(lv_path)); rv=sitk.ReadImage(str(rv_path))
    rv=sitk.Resample(rv,lv,sitk.Transform(),sitk.sitkNearestNeighbor,0,sitk.sitkUInt8)
    binary=sitk.Cast(lv>0,sitk.sitkUInt8); distance=sitk.SignedMaurerDistanceMap(binary,insideIsPositive=True,useImageSpacing=True)
    shell=(sitk.GetArrayFromImage(distance)>=0)&(sitk.GetArrayFromImage(distance)<=shell_thickness_mm)
    indices,points=_points(lv,shell); _,rvpoints=_points(lv,sitk.GetArrayFromImage(rv)>0)
    if len(points)<100 or len(rvpoints)<10: raise ValueError("LV shell or RV mask is empty/too small")
    center=points.mean(0); centered=points-center; _,vectors=np.linalg.eigh(np.cov(centered,rowvar=False)); axis=vectors[:,-1]
    # Orient long axis so the smaller cross-sectional end becomes the apex.
    projection=centered@axis
    if np.mean(np.linalg.norm(centered[projection<np.quantile(projection,.1)],axis=1)) > np.mean(np.linalg.norm(centered[projection>np.quantile(projection,.9)],axis=1)): axis=-axis
    tree=cKDTree(rvpoints); distances,_=tree.query(points); landmark=points[np.argmin(distances)]-center
    landmark-=axis*np.dot(landmark,axis); x=landmark/np.linalg.norm(landmark); y=np.cross(axis,x); y/=np.linalg.norm(y)
    xx=centered@x; yy=centered@y; zz=centered@axis; angle=np.degrees(np.arctan2(yy,xx))%360
    labels=assign_segments(zz,angle); volume=np.zeros(shell.shape,np.uint8); volume[tuple(indices.T)]=labels
    result=sitk.GetImageFromArray(volume); result.CopyInformation(lv); sitk.WriteImage(result,str(output))
    return output

