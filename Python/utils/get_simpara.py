import numpy as np

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False


def get_simpara(s_pat, sim_pat, pat_arr, para):
    f = para.f
    ch = para.ch
    l = s_pat.shape[3]
    ns = para.ns
    nsig = para.nSig

    use_cuda = _HAS_CUPY and isinstance(s_pat, cp.ndarray)
    if use_cuda:
        neighbors = s_pat[:, :, :, pat_arr]
        coe = neighbors - sim_pat[:, :, :, None, :]
        cu0 = cp.mean(coe ** 2, axis=3)
        b0 = s_pat - sim_pat
        cu0 = cp.maximum(0.0, cu0 - nsig ** 2)
        b0 = (cu0 * (b0 ** 2)) ** 0.25
        b = b0 ** 2 + para.Eps
        return b

    cu0 = np.zeros((f, f, ch, l), dtype=np.float64)
    b0 = np.zeros((f, f, ch, l), dtype=np.float64)

    for i in range(l):
        coe = s_pat[:, :, :, pat_arr[:, i]] - np.tile(sim_pat[:, :, :, i:i + 1], (1, 1, 1, ns))
        cu0[:, :, :, i] = np.mean(coe ** 2, axis=3)
        b0[:, :, :, i] = s_pat[:, :, :, i] - sim_pat[:, :, :, i]

    cu0 = np.maximum(0.0, cu0 - nsig ** 2)
    b0 = (cu0 * (b0 ** 2)) ** 0.25
    b = b0 ** 2 + para.Eps

    return b
