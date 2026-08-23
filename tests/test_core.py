import numpy as np
from heartdelta.aha import assign_segments
from heartdelta.metrics import eqd2_factor,eqd2_values,geud

def test_segment_ranges():
    z=np.r_[np.repeat(.9,6),np.repeat(.5,6),np.repeat(.2,4),0]; angle=np.r_[np.arange(0,360,60),np.arange(0,360,60),[0,90,180,270],0]
    assert set(assign_segments(z,angle))==set(range(1,18))
def test_eqd2_identity(): assert eqd2_factor(50,25,3)==1
def test_geud_mean_for_a1(): assert geud(np.array([1,2,3]),1)==2
def test_eqd2_is_converted_voxelwise_before_averaging():
    values=eqd2_values(np.array([0.0,50.0]),25,2.0)
    assert np.allclose(values,[0.0,50.0])
    assert not np.isclose(eqd2_values(np.array([0.0,25.0]),25,2.0).mean(),eqd2_values(np.array([12.5]),25,2.0).mean())
