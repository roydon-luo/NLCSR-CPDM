import numpy as np
from utils.matlab_compat import padarray
from utils.get_cpdmask import get_cpdmask
from utils.functions.bayer_residual import bayer_residual
from utils.functions.bayer_bilinear import bayer_bilinear
from utils.functions.PCDP import PCDP
from utils.functions.EARI import EARI


def init_interp(cpfa):
    h, w = cpfa.shape[0], cpfa.shape[1]
    pattern = "rggb"
    col_method = "RI"
    pol_method = "PCDP"
    mask = get_cpdmask(h, w, pattern)

    bayer_90 = cpfa[0::2, 0::2]
    bayer_45 = cpfa[0::2, 1::2]
    bayer_135 = cpfa[1::2, 0::2]
    bayer_0 = cpfa[1::2, 1::2]

    if col_method == "RI":
        sigma = 1
        eps = 1e-32
        bayer_dem_90 = bayer_residual(np.repeat(bayer_90[:, :, None], 3, axis=2), pattern, sigma, eps)
        bayer_dem_45 = bayer_residual(np.repeat(bayer_45[:, :, None], 3, axis=2), pattern, sigma, eps)
        bayer_dem_135 = bayer_residual(np.repeat(bayer_135[:, :, None], 3, axis=2), pattern, sigma, eps)
        bayer_dem_0 = bayer_residual(np.repeat(bayer_0[:, :, None], 3, axis=2), pattern, sigma, eps)
    elif col_method == "BI":
        eps = 1e-32
        mask_bayer = np.zeros((h // 2, w // 2, 3), dtype=np.float64)
        mask_bayer[0::2, 0::2, 0] = 1
        mask_bayer[0::2, 1::2, 1] = 1
        mask_bayer[1::2, 0::2, 1] = 1
        mask_bayer[1::2, 1::2, 2] = 1

        mosaic_bayer90 = bayer_90 * mask_bayer
        mosaic_bayer45 = bayer_45 * mask_bayer
        mosaic_bayer135 = bayer_135 * mask_bayer
        mosaic_bayer0 = bayer_0 * mask_bayer

        bayer_dem_90 = bayer_bilinear(mosaic_bayer90)
        bayer_dem_45 = bayer_bilinear(mosaic_bayer45)
        bayer_dem_135 = bayer_bilinear(mosaic_bayer135)
        bayer_dem_0 = bayer_bilinear(mosaic_bayer0)
    else:
        raise ValueError("Unsupported color method")

    bayer_dem_rgb = np.zeros((h, w, 3), dtype=np.float64)
    bayer_dem_rgb[0::2, 0::2, :] = bayer_dem_90
    bayer_dem_rgb[0::2, 1::2, :] = bayer_dem_45
    bayer_dem_rgb[1::2, 0::2, :] = bayer_dem_135
    bayer_dem_rgb[1::2, 1::2, :] = bayer_dem_0

    if pol_method == "PCDP":
        mask_dofp = np.zeros((h, w, 4), dtype=np.float64)
        mask_dofp[0::2, 0::2, 0] = 1
        mask_dofp[0::2, 1::2, 1] = 1
        mask_dofp[1::2, 0::2, 3] = 1
        mask_dofp[1::2, 1::2, 2] = 1
        dem_0 = np.zeros((h, w, 3), dtype=np.float64)
        dem_45 = np.zeros((h, w, 3), dtype=np.float64)
        dem_90 = np.zeros((h, w, 3), dtype=np.float64)
        dem_135 = np.zeros((h, w, 3), dtype=np.float64)
        for k in range(3):
            mosaic = padarray(bayer_dem_rgb[:, :, k][:, :, None] * mask_dofp, [8, 8, 0], mode="symmetric")
            mask_dofp_in = padarray(mask_dofp, [8, 8, 0], mode="symmetric")
            i90, i45, i0, i135 = PCDP(mosaic, mask_dofp_in)
            dem_0[:, :, k] = i0
            dem_45[:, :, k] = i45
            dem_90[:, :, k] = i90
            dem_135[:, :, k] = i135
    elif pol_method == "EARI":
        mask_p90 = mask[:, :, 6] + mask[:, :, 7] + mask[:, :, 8]
        mask_p45 = mask[:, :, 3] + mask[:, :, 4] + mask[:, :, 5]
        mask_p135 = mask[:, :, 9] + mask[:, :, 10] + mask[:, :, 11]
        mask_p0 = mask[:, :, 0] + mask[:, :, 1] + mask[:, :, 2]
        dem_0, dem_45, dem_90, dem_135 = EARI(bayer_dem_rgb, eps, mask_p0, mask_p45, mask_p90, mask_p135)
    else:
        raise ValueError("Unsupported polarization method")

    s_ini = np.concatenate([dem_0, dem_45, dem_90, dem_135], axis=2)
    return s_ini, mask
