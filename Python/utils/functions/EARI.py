import numpy as np
from utils.matlab_compat import imfilter
from utils.functions.residual_interpolation import residual_interpolation


def EARI(mpfa, eps, mask_p0, mask_p45, mask_p90, mask_p135):
    polar_90 = np.zeros_like(mpfa, dtype=np.float64)
    polar_45 = np.zeros_like(mpfa, dtype=np.float64)
    polar_135 = np.zeros_like(mpfa, dtype=np.float64)
    polar_0 = np.zeros_like(mpfa, dtype=np.float64)

    fn = np.array([
        [1 / 8, 1 / 4, 1 / 8],
        [1 / 8, 1 / 4, 1 / 8],
        [0, 0, 0],
    ], dtype=np.float64)
    fs = np.array([
        [0, 0, 0],
        [1 / 8, 1 / 4, 1 / 8],
        [1 / 8, 1 / 4, 1 / 8],
    ], dtype=np.float64)
    fw = fn.T
    fe = fs.T

    xn = imfilter(mpfa, fn, "replicate")
    xe = imfilter(mpfa, fe, "replicate")
    xw = imfilter(mpfa, fw, "replicate")
    xs = imfilter(mpfa, fs, "replicate")

    hn = np.array([
        [-1 / 2, 1, -1 / 2],
        [1 / 2, -1, 1 / 2],
        [0, 0, 0],
    ], dtype=np.float64)
    hs = np.array([
        [0, 0, 0],
        [1 / 2, -1, 1 / 2],
        [-1 / 2, 1, -1 / 2],
    ], dtype=np.float64)
    hw = hn.T
    he = hs.T

    inn = np.abs(imfilter(mpfa, hn, "replicate"))
    ie = np.abs(imfilter(mpfa, he, "replicate"))
    iw = np.abs(imfilter(mpfa, hw, "replicate"))
    is_ = np.abs(imfilter(mpfa, hs, "replicate"))

    mn = np.array([
        [1 / 15] * 5,
        [1 / 15] * 5,
        [1 / 15] * 5,
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], dtype=np.float64)
    ms = np.array([
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1 / 15] * 5,
        [1 / 15] * 5,
        [1 / 15] * 5,
    ], dtype=np.float64)
    mw = mn.T
    me = ms.T

    wn = imfilter(inn, mn, "replicate")
    we = imfilter(ie, me, "replicate")
    ww = imfilter(iw, mw, "replicate")
    ws = imfilter(is_, ms, "replicate")

    wn = 1.0 / (wn + eps)
    we = 1.0 / (we + eps)
    ww = 1.0 / (ww + eps)
    ws = 1.0 / (ws + eps)

    wsum = wn + we + ww + ws
    guide = (wn * xn + we * xe + ww * xw + ws * xs) / (wsum + eps)

    for c in range(mpfa.shape[2]):
        polar_90[:, :, c], _ = residual_interpolation(guide[:, :, c], mask_p90 * mpfa[:, :, c], mask_p90, eps)
        polar_45[:, :, c], _ = residual_interpolation(guide[:, :, c], mask_p45 * mpfa[:, :, c], mask_p45, eps)
        polar_135[:, :, c], _ = residual_interpolation(guide[:, :, c], mask_p135 * mpfa[:, :, c], mask_p135, eps)
        polar_0[:, :, c], _ = residual_interpolation(guide[:, :, c], mask_p0 * mpfa[:, :, c], mask_p0, eps)

    return polar_0, polar_45, polar_90, polar_135
