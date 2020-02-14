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

    assert np.all((mis - np.array([
        0.16816107, 0.16816107, 0.02576454, 0.16816107, 0.16816107, 0.02576454,
        0.16816107, 0.16816107, 0.02576454, 0.42107727, 0.42107727, 0.42107727
    ]) < 1e-7))

    hitsA, hitsB, mis = spydrpick(pa_matrix,
                                  weights=[0, 1, 1, 1, 1, 1],
                                  keep_quantile=0,
                                  chunk_size=100)

    assert np.all((mis - np.array([
        0.14913948, 0.14913948, -0.01359942, 0.14913948, 0.14913948,
        -0.01359942, 0.14913948, 0.14913948, -0.01359942, 0.24591178,
        0.24591178, 0.24591178
    ]) < 1e-7))

    return
