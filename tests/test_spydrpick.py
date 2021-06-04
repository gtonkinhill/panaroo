# test the central spydrpick function
from panaroo.spydrpick import spydrpick
import numpy as np


def test_spydrpick(datafolder):

    pa_matrix = np.array([[1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0],
                            [0, 0, 0, 1, 1, 1], [1, 1, 1, 1, 1, 0]])

    hitsA, hitsB, mis = zip(*spydrpick(pa_matrix,
                                    weights=None,
                                    keep_quantile=0,
                                    chunk_size=100))

    assert np.all((np.array(mis) - np.array([
        0.316377019, 0.04316844491, 0.04316844491, 
        0.316377019, 0.04316844491, 0.3163770193
    ]) < 1e-7))

    return
