import numpy as np
from utils.matlab_compat import conv2


def PCDP(mosaic, mask):
    weight_orth = 1.0 / (1.0 + 2.0 * np.sqrt(2.0))
    weight_no_orth = np.sqrt(2.0) / (1.0 + 2.0 * np.sqrt(2.0))

    conv_kernel = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], dtype=np.float64) / 4.0

    i0 = conv2(mosaic[:, :, 0], conv_kernel, mode="same")
    i45 = conv2(mosaic[:, :, 1], conv_kernel, mode="same")
    i90 = conv2(mosaic[:, :, 2], conv_kernel, mode="same")
    i135 = conv2(mosaic[:, :, 3], conv_kernel, mode="same")

    bi = np.stack([i0, i45, i90, i135], axis=2)
    array = [4, 1, 2, 3, 4, 1, 2]
    demosaic = np.zeros((mosaic.shape[0], mosaic.shape[1], 4), dtype=np.float64)

    for k in range(4):
        idx = array[k]
        i_diff_1 = mosaic[:, :, array[k + 1] - 1] - bi[:, :, idx - 1] * mask[:, :, array[k + 1] - 1]
        i_diff_2 = mosaic[:, :, array[k + 1] - 1] - bi[:, :, array[k + 2] - 1] * mask[:, :, array[k + 1] - 1]
        i_diff_3 = mosaic[:, :, array[k + 1] - 1] - bi[:, :, array[k + 3] - 1] * mask[:, :, array[k + 1] - 1]

        i_1 = conv2(i_diff_1, conv_kernel, mode="same")
        i_2 = conv2(i_diff_2, conv_kernel, mode="same")
        i_3 = conv2(i_diff_3, conv_kernel, mode="same")

        i_1 = bi[:, :, idx - 1] + i_1
        i_2 = bi[:, :, array[k + 2] - 1] + i_2
        i_3 = bi[:, :, array[k + 3] - 1] + i_3

        demosaic[:, :, array[k + 1] - 1] = weight_no_orth * (i_1 + i_2) + weight_orth * i_3

    demosaic = demosaic[8:-8, 8:-8, :]
    i0 = demosaic[:, :, 0]
    i45 = demosaic[:, :, 1]
    i90 = demosaic[:, :, 2]
    i135 = demosaic[:, :, 3]

    # Match MATLAB output order [i90, i45, i0, i135]
    i90_out = i0
    i0_out = i90
    return i90_out, i45, i0_out, i135
