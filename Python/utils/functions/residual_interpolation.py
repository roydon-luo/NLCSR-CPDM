import numpy as np
from utils.matlab_compat import imfilter
from utils.functions.guidedfilter_MLRI_wei import guidedfilter_MLRI_wei
from utils.functions.clip import clip


def residual_interpolation(guide, mosaic, mask, eps):
    h = 5
    v = 5

    f = np.array(
        [
            [0, 0, -1, 0, 0],
            [0, 0, 0, 0, 0],
            [-1, 0, 4, 0, -1],
            [0, 0, 0, 0, 0],
            [0, 0, -1, 0, 0],
        ],
        dtype=np.float64,
    )
    lap_input = imfilter(mosaic, f, "replicate")
    lap_guide = imfilter(guide * mask, f, "replicate")

    tentative, dif = guidedfilter_MLRI_wei(guide, mosaic, mask, lap_guide, lap_input, mask, h, v, eps)
    tentative = clip(tentative, 0, 1)
    residual = mask * (mosaic - tentative)

    hker = np.array(
        [
            [1 / 4, 1 / 2, 1 / 4],
            [1 / 2, 1, 1 / 2],
            [1 / 4, 1 / 2, 1 / 4],
        ],
        dtype=np.float64,
    )
    residual = imfilter(residual, hker, "replicate")
    output = residual + tentative

    return output, dif
