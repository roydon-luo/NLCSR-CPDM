import numpy as np
from scipy.spatial.distance import cdist

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False


def get_simblk(s_pat, para):
    f = para.f
    ch = para.ch
    l = s_pat.shape[3]
    ns = para.ns

    use_cuda = _HAS_CUPY and isinstance(s_pat, cp.ndarray)
    if use_cuda:
        s_vec = cp.reshape(s_pat, (f * f * ch, l), order="F").astype(cp.float64, copy=False)
        norms = cp.sum(s_vec * s_vec, axis=0)
        gram = s_vec.T @ s_vec
        dist2 = norms[:, None] + norms[None, :] - 2.0 * gram
        dist2 = cp.maximum(dist2, 0.0)
        distances = cp.sqrt(dist2)
        distances[cp.arange(l), cp.arange(l)] = cp.inf

        sorted_ind = cp.argsort(distances, axis=1)
        sorted_dis = cp.take_along_axis(distances, sorted_ind, axis=1)
        min_dis = sorted_dis[:, 0:ns]
        min_ind = sorted_ind[:, 0:ns]

        pat_arr = min_ind.T
        valid_min_dis = (min_dis < para.min_dis).astype(cp.float64)
        valid_dis = valid_min_dis * min_dis
        wei = cp.exp(-valid_dis / 10.0)
        wei_arr = wei / (cp.sum(wei, axis=1, keepdims=True) + cp.finfo(cp.float64).eps)
        wei_arr = wei_arr.T

        sim_pat = cp.zeros((f, f, ch, l), dtype=cp.float64)
        for i in range(ns):
            w = wei_arr[i, :].reshape(1, 1, 1, l)
            sim_pat = sim_pat + s_pat[:, :, :, pat_arr[i, :]] * w
        return sim_pat, pat_arr

    s_vec = np.reshape(s_pat, (f * f * ch, l), order="F")
    distances = cdist(s_vec.T, s_vec.T, metric="euclidean")
    np.fill_diagonal(distances, np.inf)

    # MATLAB sort is stable for ties; keep neighbor order deterministic.
    sorted_ind = np.argsort(distances, axis=1, kind="stable")
    sorted_dis = np.take_along_axis(distances, sorted_ind, axis=1)
    min_dis = sorted_dis[:, 0:ns]
    min_ind = sorted_ind[:, 0:ns]

    pat_arr = min_ind.T
    valid_min_dis = (min_dis < para.min_dis).astype(np.float64)
    valid_dis = valid_min_dis * min_dis
    wei = np.exp(-valid_dis / 10.0)
    wei_arr = wei / (np.sum(wei, axis=1, keepdims=True) + np.finfo(np.float64).eps)
    wei_arr = wei_arr.T

    sim_pat = np.zeros((f, f, ch, l), dtype=np.float64)
    for i in range(ns):
        w = wei_arr[i, :].reshape(1, 1, 1, l)
        w = np.tile(w, (f, f, ch, 1))
        sim_cur = s_pat[:, :, :, pat_arr[i, :]] * w
        sim_pat = sim_pat + sim_cur

    return sim_pat, pat_arr
