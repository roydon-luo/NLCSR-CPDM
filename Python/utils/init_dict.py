import numpy as np
from utils.matlab_compat import padarray


def init_dict(para):
    m = np.atleast_1d(para.m)
    k = np.atleast_1d(para.K)
    s = np.max(m)
    d_all = []

    for i in range(len(m)):
        d = np.random.randn(int(m[i]), int(m[i]), int(k[i])).astype(np.float64)
        denom = np.sqrt(np.sum(d ** 2, axis=(0, 1), keepdims=True))
        d = d / denom
        d = padarray(d, [int(s - m[i]), int(s - m[i])], mode="post")
        d_all.append(d)

    if len(d_all) == 1:
        return d_all[0]
    return np.stack(d_all, axis=3)
