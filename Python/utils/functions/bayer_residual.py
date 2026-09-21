import numpy as np
from utils.functions.mosaic_bayer import mosaic_bayer
from utils.functions.green_interpolation import green_interpolation
from utils.functions.red_interpolation import red_interpolation
from utils.functions.blue_interpolation import blue_interpolation
from utils.functions.clip import clip


def bayer_residual(rgb, pattern, sigma, eps):
    mosaic, mask = mosaic_bayer(rgb, pattern)

    green = green_interpolation(mosaic, mask, pattern, sigma, eps)
    green = clip(green, 0, 255)

    red = red_interpolation(green, mosaic, mask, eps)
    blue = blue_interpolation(green, mosaic, mask, eps)
    red = clip(red, 0, 255)
    blue = clip(blue, 0, 255)

    rgb_dem = np.zeros_like(rgb, dtype=np.float64)
    rgb_dem[:, :, 0] = red
    rgb_dem[:, :, 1] = green
    rgb_dem[:, :, 2] = blue

    return rgb_dem
