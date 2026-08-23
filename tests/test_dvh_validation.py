import numpy as np
from heartdelta.dvh import structure_metrics
from heartdelta.validation import dice


def test_dvh_metrics():
    out=structure_metrics(np.array([0.,10.,20.,30.]),.01)
    assert out["mean_dose_eqd2"]==15
    assert out["v20_percent"]==50


def test_dice_identity():
    mask=np.array([0,1,1,0])
    assert dice(mask,mask)==1
