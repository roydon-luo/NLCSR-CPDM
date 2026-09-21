import numpy as np
from utils.matlab_compat import imfilter, fspecial_gaussian
from utils.functions.guidedfilter_MLRI_wei import guidedfilter_MLRI_wei
from utils.functions.clip import clip


def green_interpolation(mosaic, mask, pattern, sigma, eps):
    imask = (mask == 0).astype(np.float64)
    rawq = np.sum(mosaic, axis=2)

    mask_gr = np.zeros_like(rawq)
    mask_gb = np.zeros_like(rawq)
    if pattern == "grbg":
        mask_gr[0::2, 0::2] = 1
        mask_gb[1::2, 1::2] = 1
    elif pattern == "rggb":
        mask_gr[0::2, 1::2] = 1
        mask_gb[1::2, 0::2] = 1
    elif pattern == "gbrg":
        mask_gb[0::2, 0::2] = 1
        mask_gr[1::2, 1::2] = 1
    elif pattern == "bggr":
        mask_gb[0::2, 1::2] = 1
        mask_gr[1::2, 0::2] = 1

    kh = np.array([0.5, 0.0, 0.5], dtype=np.float64)
    kv = kh.reshape(-1, 1)
    rawh = imfilter(rawq, kh, "replicate")
    rawv = imfilter(rawq, kv, "replicate")
    guide_gh = mosaic[:, :, 1] + rawh * mask[:, :, 0] + rawh * mask[:, :, 2]
    guide_rh = mosaic[:, :, 0] + rawh * mask_gr
    guide_bh = mosaic[:, :, 2] + rawh * mask_gb
    guide_gv = mosaic[:, :, 1] + rawv * mask[:, :, 0] + rawv * mask[:, :, 2]
    guide_rv = mosaic[:, :, 0] + rawv * mask_gb
    guide_bv = mosaic[:, :, 2] + rawv * mask_gr

    h = 3
    v = 3
    f = np.array([-1, 0, 2, 0, -1], dtype=np.float64)

    dif_r = imfilter(mosaic[:, :, 0], f, "replicate")
    dif_gr = imfilter(guide_gh * mask[:, :, 0], f, "replicate")
    tentative_rh, _ = guidedfilter_MLRI_wei(guide_gh, mosaic[:, :, 0], mask[:, :, 0], dif_gr, dif_r, mask[:, :, 0], h, v, eps)

    dif_gr = imfilter(mosaic[:, :, 1] * mask_gr, f, "replicate")
    dif_r = imfilter(guide_rh * mask_gr, f, "replicate")
    tentative_grh, _ = guidedfilter_MLRI_wei(guide_rh, mosaic[:, :, 1] * mask_gr, mask_gr, dif_r, dif_gr, mask_gr, h, v, eps)

    dif_b = imfilter(mosaic[:, :, 2], f, "replicate")
    dif_gb = imfilter(guide_gh * mask[:, :, 2], f, "replicate")
    tentative_bh, _ = guidedfilter_MLRI_wei(guide_gh, mosaic[:, :, 2], mask[:, :, 2], dif_gb, dif_b, mask[:, :, 2], h, v, eps)

    dif_gb = imfilter(mosaic[:, :, 1] * mask_gb, f, "replicate")
    dif_b = imfilter(guide_bh * mask_gb, f, "replicate")
    tentative_gbh, _ = guidedfilter_MLRI_wei(guide_bh, mosaic[:, :, 1] * mask_gb, mask_gb, dif_b, dif_gb, mask_gb, h, v, eps)

    f = f.reshape(-1, 1)
    dif_r = imfilter(mosaic[:, :, 0], f, "replicate")
    dif_gr = imfilter(guide_gv * mask[:, :, 0], f, "replicate")
    tentative_rv, _ = guidedfilter_MLRI_wei(guide_gv, mosaic[:, :, 0], mask[:, :, 0], dif_gr, dif_r, mask[:, :, 0], v, h, eps)

    dif_gr = imfilter(mosaic[:, :, 1] * mask_gb, f, "replicate")
    dif_r = imfilter(guide_rv * mask_gb, f, "replicate")
    tentative_grv, _ = guidedfilter_MLRI_wei(guide_rv, mosaic[:, :, 1] * mask_gb, mask_gb, dif_r, dif_gr, mask_gb, v, h, eps)

    dif_b = imfilter(mosaic[:, :, 2], f, "replicate")
    dif_gb = imfilter(guide_gv * mask[:, :, 2], f, "replicate")
    tentative_bv, _ = guidedfilter_MLRI_wei(guide_gv, mosaic[:, :, 2], mask[:, :, 2], dif_gb, dif_b, mask[:, :, 2], v, h, eps)

    dif_gb = imfilter(mosaic[:, :, 1] * mask_gr, f, "replicate")
    dif_b = imfilter(guide_bv * mask_gr, f, "replicate")
    tentative_gbv, _ = guidedfilter_MLRI_wei(guide_bv, mosaic[:, :, 1] * mask_gr, mask_gr, dif_b, dif_gb, mask_gr, v, h, eps)

    tentative_grh = clip(tentative_grh, 0, 255)
    tentative_grv = clip(tentative_grv, 0, 255)
    tentative_gbh = clip(tentative_gbh, 0, 255)
    tentative_gbv = clip(tentative_gbv, 0, 255)
    tentative_rh = clip(tentative_rh, 0, 255)
    tentative_rv = clip(tentative_rv, 0, 255)
    tentative_bh = clip(tentative_bh, 0, 255)
    tentative_bv = clip(tentative_bv, 0, 255)

    residual_grh = (mosaic[:, :, 1] - tentative_grh) * mask_gr
    residual_gbh = (mosaic[:, :, 1] - tentative_gbh) * mask_gb
    residual_rh = (mosaic[:, :, 0] - tentative_rh) * mask[:, :, 0]
    residual_bh = (mosaic[:, :, 2] - tentative_bh) * mask[:, :, 2]
    residual_grv = (mosaic[:, :, 1] - tentative_grv) * mask_gb
    residual_gbv = (mosaic[:, :, 1] - tentative_gbv) * mask_gr
    residual_rv = (mosaic[:, :, 0] - tentative_rv) * mask[:, :, 0]
    residual_bv = (mosaic[:, :, 2] - tentative_bv) * mask[:, :, 2]

    kh = np.array([0.5, 0.0, 0.5], dtype=np.float64)
    residual_grh = imfilter(residual_grh, kh, "replicate")
    residual_gbh = imfilter(residual_gbh, kh, "replicate")
    residual_rh = imfilter(residual_rh, kh, "replicate")
    residual_bh = imfilter(residual_bh, kh, "replicate")
    kv = kh.reshape(-1, 1)
    residual_grv = imfilter(residual_grv, kv, "replicate")
    residual_gbv = imfilter(residual_gbv, kv, "replicate")
    residual_rv = imfilter(residual_rv, kv, "replicate")
    residual_bv = imfilter(residual_bv, kv, "replicate")

    grh = (tentative_grh + residual_grh) * mask[:, :, 0]
    gbh = (tentative_gbh + residual_gbh) * mask[:, :, 2]
    rh = (tentative_rh + residual_rh) * mask_gr
    bh = (tentative_bh + residual_bh) * mask_gb
    grv = (tentative_grv + residual_grv) * mask[:, :, 0]
    gbv = (tentative_gbv + residual_gbv) * mask[:, :, 2]
    rv = (tentative_rv + residual_rv) * mask_gb
    bv = (tentative_bv + residual_bv) * mask_gr

    grh = clip(grh, 0, 255)
    grv = clip(grv, 0, 255)
    gbh = clip(gbh, 0, 255)
    gbv = clip(gbv, 0, 255)
    rh = clip(rh, 0, 255)
    rv = clip(rv, 0, 255)
    bh = clip(bh, 0, 255)
    bv = clip(bv, 0, 255)

    difh = mosaic[:, :, 1] + grh + gbh - mosaic[:, :, 0] - mosaic[:, :, 2] - rh - bh
    difv = mosaic[:, :, 1] + grv + gbv - mosaic[:, :, 0] - mosaic[:, :, 2] - rv - bv

    kh = np.array([1, 0, -1], dtype=np.float64)
    kv = kh.reshape(-1, 1)
    difh2 = np.abs(imfilter(difh, kh, "replicate"))
    difv2 = np.abs(imfilter(difv, kv, "replicate"))

    k = np.ones((3, 3), dtype=np.float64)
    wh = imfilter(difh2, k, "replicate")
    wv = imfilter(difv2, k, "replicate")
    kw = np.array([1, 0, 0], dtype=np.float64)
    ke = np.array([0, 0, 1], dtype=np.float64)
    ks = ke.reshape(-1, 1)
    kn = kw.reshape(-1, 1)
    ww = imfilter(wh, kw, "replicate")
    we = imfilter(wh, ke, "replicate")
    wn = imfilter(wv, kn, "replicate")
    ws = imfilter(wv, ks, "replicate")
    ww = 1.0 / (np.power(ww, 2) + 1e-2)
    we = 1.0 / (np.power(we, 2) + 1e-2)
    ws = 1.0 / (np.power(ws, 2) + 1e-2)
    wn = 1.0 / (np.power(wn, 2) + 1e-2)

    hker = fspecial_gaussian([1, 9], sigma)
    ke = np.array([0, 0, 0, 0, 1, 1, 1, 1, 1], dtype=np.float64) * hker
    kw = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=np.float64) * hker
    ke = ke / np.sum(ke)
    kw = kw / np.sum(kw)
    ks = ke.reshape(-1, 1)
    kn = kw.reshape(-1, 1)
    difn = imfilter(difv, kn, "replicate")
    difs = imfilter(difv, ks, "replicate")
    difw = imfilter(difh, kw, "replicate")
    dife = imfilter(difh, ke, "replicate")
    wt = ww + we + wn + ws
    dif = (wn * difn + ws * difs + ww * difw + we * dife) / wt
    green = dif + rawq
    green = green * imask[:, :, 1] + rawq * mask[:, :, 1]

    return green
