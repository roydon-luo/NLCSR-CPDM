import numpy as np
from utils.functions.boxfilter import boxfilter


def guidedfilter_MLRI_wei(g, r, mask, i, p, m, h, v, eps):
    hei, wid = i.shape
    n = boxfilter(m, h, v)
    n[n == 0] = 1
    _ = boxfilter(np.ones((hei, wid)), h, v)

    mean_ip = boxfilter(i * p * m, h, v) / n
    mean_ii = boxfilter(i * i * m, h, v) / n

    a = mean_ip / (mean_ii + eps)

    n3 = boxfilter(mask, h, v)
    n3[n3 == 0] = 1
    mean_g = boxfilter(g * mask, h, v) / n3
    mean_r = boxfilter(r * mask, h, v) / n3
    b = mean_r - a * mean_g

    dif2 = (
        boxfilter(g * g * mask, h, v) * a * a
        + b * b * n3
        + boxfilter(r * r * mask, h, v)
        + 2 * a * b * boxfilter(g * mask, h, v)
        - 2 * b * boxfilter(r * mask, h, v)
        - 2 * a * boxfilter(r * g * mask, h, v)
    )
    dif2 = dif2 / n3
    dif = dif2.copy()
    dif[dif < 0.01] = 0.01
    dif = 1.0 / dif
    wdif = boxfilter(dif, h, v)
    wdif[wdif < 0.01] = 0.01
    mean_a = boxfilter(a * dif, h, v) / wdif
    mean_b = boxfilter(b * dif, h, v) / wdif

    q = mean_a * g + mean_b
    return q, dif2
