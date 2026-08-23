import numpy as np
from heartdelta.aha import assign_segments
from heartdelta.metrics import eqd2_factor,geud

def test_segment_ranges():
    z=np.r_[np.repeat(.9,6),np.repeat(.5,6),np.repeat(.2,4),0]; angle=np.r_[np.arange(0,360,60),np.arange(0,360,60),[0,90,180,270],0]
    assert set(assign_segments(z,angle))==set(range(1,18))
def test_eqd2_identity(): assert eqd2_factor(50,25,3)==1
def test_geud_mean_for_a1(): assert geud(np.array([1,2,3]),1)==2

