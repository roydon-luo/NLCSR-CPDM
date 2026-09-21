import numpy as np
from utils.matlab_compat import imfilter
from utils.functions.guidedfilter_MLRI_wei import guidedfilter_MLRI_wei
from utils.functions.clip import clip


def blue_interpolation(green, mosaic, mask, eps):
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
    lap_blue = imfilter(mosaic[:, :, 2], f, "replicate")
    lap_green = imfilter(green * mask[:, :, 2], f, "replicate")

    tentative_b, _ = guidedfilter_MLRI_wei(
        green, mosaic[:, :, 2], mask[:, :, 2], lap_green, lap_blue, mask[:, :, 2], h, v, eps
    )
    tentative_b = clip(tentative_b, 0, 255)
    residual_b = mask[:, :, 2] * (mosaic[:, :, 2] - tentative_b)

    hker = np.array(
        [
            [1 / 4, 1 / 2, 1 / 4],
            [1 / 2, 1, 1 / 2],
            [1 / 4, 1 / 2, 1 / 4],
        ],
        dtype=np.float64,
    )
    residual_b = imfilter(residual_b, hker, "replicate")
    blue = residual_b + tentative_b

    return blue
