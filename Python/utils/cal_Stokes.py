import numpy as np


def cal_Stokes(im):
    im_out0 = im[:, :, 0:3]
    im_out45 = im[:, :, 3:6]
    im_out90 = im[:, :, 6:9]
    im_out135 = im[:, :, 9:12]

    s0 = (im_out0 + im_out45 + im_out90 + im_out135) * 0.5
    s1 = im_out0 - im_out90
    s2 = im_out45 - im_out135
    dolp = np.sqrt(s1 ** 2 + s2 ** 2) / (s0 + np.finfo(np.float64).eps)
    aolp = 0.5 * np.arctan2(s2, s1)

    return s0, dolp, aolp
