import numpy as np
from utils.matlab_compat import conv2


def bayer_bilinear(mosaic):
    conv_kernel_g = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], dtype=np.float64) / 8.0
    conv_kernel_r = np.array([[1, 0, 1], [0, 4, 0], [1, 0, 1]], dtype=np.float64) / 4.0

    r_1 = conv2(mosaic[:, :, 0], conv_kernel_r, mode="same")
    r = conv2(r_1, conv_kernel_g, mode="same")
    g = conv2(mosaic[:, :, 1], conv_kernel_g, mode="same")
    b_1 = conv2(mosaic[:, :, 2], conv_kernel_r, mode="same")
    b = conv2(b_1, conv_kernel_g, mode="same")

    return np.stack([r, g, b], axis=2)
