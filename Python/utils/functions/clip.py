import numpy as np


def clip(x, lo, hi):
    y = np.array(x, dtype=np.float64, copy=True)
    y[y < lo] = lo
    y[y > hi] = hi
    return y
