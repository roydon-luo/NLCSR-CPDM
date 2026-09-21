import numpy as np


def mosaic_bayer(rgb, pattern):
    num = np.zeros(len(pattern), dtype=np.int64)
    p = [i for i, ch in enumerate(pattern) if ch in ("r", "R")]
    for idx in p:
        num[idx] = 1
    p = [i for i, ch in enumerate(pattern) if ch in ("g", "G")]
    for idx in p:
        num[idx] = 2
    p = [i for i, ch in enumerate(pattern) if ch in ("b", "B")]
    for idx in p:
        num[idx] = 3

    h, w = rgb.shape[0], rgb.shape[1]
    mosaic = np.zeros((h, w, 3), dtype=np.float64)
    mask = np.zeros((h, w, 3), dtype=np.float64)

    rows1 = np.arange(0, h, 2)
    rows2 = np.arange(1, h, 2)
    cols1 = np.arange(0, w, 2)
    cols2 = np.arange(1, w, 2)

    mask[np.ix_(rows1, cols1, [num[0] - 1])] = 1
    mask[np.ix_(rows1, cols2, [num[1] - 1])] = 1
    mask[np.ix_(rows2, cols1, [num[2] - 1])] = 1
    mask[np.ix_(rows2, cols2, [num[3] - 1])] = 1

    mosaic[:, :, 0] = rgb[:, :, 0] * mask[:, :, 0]
    mosaic[:, :, 1] = rgb[:, :, 1] * mask[:, :, 1]
    mosaic[:, :, 2] = rgb[:, :, 2] * mask[:, :, 2]

    return mosaic, mask
