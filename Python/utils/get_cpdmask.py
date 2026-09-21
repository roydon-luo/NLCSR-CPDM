import numpy as np


def get_cpdmask(h, w, pattern):
    num = np.zeros(len(pattern), dtype=np.int64)
    pr = [i for i, ch in enumerate(pattern) if ch in ("r", "R")]
    for idx in pr:
        num[idx] = 1
    pg = [i for i, ch in enumerate(pattern) if ch in ("g", "G")]
    for idx in pg:
        num[idx] = 2
    pb = [i for i, ch in enumerate(pattern) if ch in ("b", "B")]
    for idx in pb:
        num[idx] = 3

    mask_rgb = np.zeros((h // 2, w // 2, 3), dtype=np.float64)
    mask_rgb[0::2, 0::2, num[0] - 1] = 1
    mask_rgb[0::2, 1::2, num[1] - 1] = 1
    mask_rgb[1::2, 0::2, num[2] - 1] = 1
    mask_rgb[1::2, 1::2, num[3] - 1] = 1

    mask = np.zeros((h, w, 12), dtype=np.float64)
    for i in range(3):
        mask[1::2, 1::2, i] = mask_rgb[:, :, i]
        mask[0::2, 1::2, i + 3] = mask_rgb[:, :, i]
        mask[0::2, 0::2, i + 6] = mask_rgb[:, :, i]
        mask[1::2, 0::2, i + 9] = mask_rgb[:, :, i]

    return mask
