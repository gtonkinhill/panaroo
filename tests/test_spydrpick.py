# test the central spydrpick function
from panaroo.spydrpick import spydrpick
import numpy as np


def test_spydrpick(datafolder):

    pa_matrix = np.array([[1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0],
                            [0, 0, 0, 1, 1, 1], [1, 1, 1, 1, 1, 0]])

    hitsA, hitsB, mis = spydrpick(pa_matrix,
                                    weights=None,
                                    keep_quantile=0,
                                    chunk_size=100)

    res = np.zeros((4,4))
    res[hitsA,hitsB] = mis

    assert np.all((res - np.array([[0.        , 0.31637702, 0.31637702, 0.04316844],
       [0.        , 0.        , 0.31637702, 0.04316844],
       [0.        , 0.        , 0.        , 0.04316844],
       [0.        , 0.        , 0.        , 0.        ]]) < 1e-4))

    return
