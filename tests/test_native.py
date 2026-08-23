import numpy as np
import pandas as pd
from heartdelta.native import aggregate_levels


def test_level_aggregation_is_voxel_weighted():
    frame = pd.DataFrame([
        {"case_id":"a","timepoint":"fu1","level":"basal","segment":1,"n_voxels":1,"delta_hu_mean":10.0},
        {"case_id":"a","timepoint":"fu1","level":"basal","segment":2,"n_voxels":3,"delta_hu_mean":2.0},
    ])
    out = aggregate_levels(frame)
    assert out.loc[0,"n_voxels"] == 4
    assert np.isclose(out.loc[0,"delta_hu_mean"], 4.0)
