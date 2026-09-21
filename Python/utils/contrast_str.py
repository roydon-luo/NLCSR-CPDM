import numpy as np


def contrast_str(x, a):
    y = np.array(x, dtype=np.float64, copy=True)
    y[y > a] = a
    y[y < 0] = 0
    y = y / a
    return y
