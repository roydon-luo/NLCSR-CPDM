import numpy as np
from utils.matlab_compat import imfilter
from utils.functions.guidedfilter_MLRI_wei import guidedfilter_MLRI_wei
from utils.functions.clip import clip


def red_interpolation(green, mosaic, mask, eps):
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
    lap_red = imfilter(mosaic[:, :, 0], f, "replicate")
    lap_green = imfilter(green * mask[:, :, 0], f, "replicate")

    tentative_r, _ = guidedfilter_MLRI_wei(
        green, mosaic[:, :, 0], mask[:, :, 0], lap_green, lap_red, mask[:, :, 0], h, v, eps
    )
    tentative_r = clip(tentative_r, 0, 255)
    residual_r = mask[:, :, 0] * (mosaic[:, :, 0] - tentative_r)

    hker = np.array(
        [
            [1 / 4, 1 / 2, 1 / 4],
            [1 / 2, 1, 1 / 2],
            [1 / 4, 1 / 2, 1 / 4],
        ],
        dtype=np.float64,
    )
    residual_r = imfilter(residual_r, hker, "replicate")
    red = residual_r + tentative_r

    return red
