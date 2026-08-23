import numpy as np
import SimpleITK as sitk
from heartdelta.calibration import apply, parameters


def test_trachea_aorta_two_point_mapping(tmp_path):
    array=np.zeros((5,5,5),np.float32); array[:2]=-950; array[2:4]=40
    ct=sitk.GetImageFromArray(array); trachea=sitk.GetImageFromArray((array==-950).astype(np.uint8)); aorta=sitk.GetImageFromArray((array==40).astype(np.uint8))
    trachea.CopyInformation(ct); aorta.CopyInformation(ct); tp=tmp_path/"trachea.nii.gz"; ap=tmp_path/"aorta.nii.gz"; sitk.WriteImage(trachea,str(tp)); sitk.WriteImage(aorta,str(ap))
    result=parameters(ct,tp,ap,min_voxels=10); mapped=apply(np.array([-950.,40.]),result)
    assert np.allclose(mapped,[-1000.,50.])
