import numpy as np
from heartdelta.slab import _center, _trimmed


def test_trimmed_mean_removes_balanced_tails():
    assert _trimmed(np.arange(10.0)) == 4.5


def test_slab_center_is_stable_for_symmetric_projection():
    value = _center(np.repeat(np.arange(-4.0, 5.0), [1, 2, 3, 4, 8, 4, 3, 2, 1]))
    assert abs(value) <= 0.5
