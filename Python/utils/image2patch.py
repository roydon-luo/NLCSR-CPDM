import numpy as np

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False

def image2patch(s, para):
    use_cuda = _HAS_CUPY and isinstance(s, cp.ndarray)
    xp = cp if use_cuda else np

    h, w, ch = s.shape
    f = para.f
    step = para.step
    h1 = h - f + 1
    w1 = w - f + 1
    r = list(range(0, h1, step))
    if r[-1] != h1 - 1:
        r.extend(list(range(r[-1] + 1, h1)))
    c = list(range(0, w1, step))
    if c[-1] != w1 - 1:
        c.extend(list(range(c[-1] + 1, w1)))
    para.r = np.asarray(r, dtype=np.int32)
    para.c = np.asarray(c, dtype=np.int32)
    nr = len(r)
    nc = len(c)

    s = xp.asarray(s, dtype=xp.float64)
    s0, s1, s2 = s.strides
    win = xp.lib.stride_tricks.as_strided(
        s,
        shape=(h - f + 1, w - f + 1, f, f, ch),
        strides=(s0, s1, s0, s1, s2),
    )

    rr = xp.asarray(para.r)
    cc = xp.asarray(para.c)
    blk = win[rr[:, None], cc[None, :], :, :, :]
    s_blk = xp.transpose(blk, (2, 3, 4, 0, 1)).reshape((f, f, ch, nr * nc), order="F")

    return s_blk, para
