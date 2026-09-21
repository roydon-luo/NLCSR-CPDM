import numpy as np
import os
from utils.fft_backend import fft2 as fft2_backend, ifft2 as ifft2_backend, fft2_gpu, ifft2_gpu

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False


def _fft2(x):
    return fft2_backend(x, axes=(0, 1))


def _ifft2(x):
    return np.real(ifft2_backend(x, axes=(0, 1)))


def _use_cuda():
    return _HAS_CUPY and (os.environ.get("USE_CUDA", "0").strip() == "1")


def G_update(x, w, s, sim, sig, eta, para, return_debug=False, s_f=None, sim_f=None):
    use_cuda = _HAS_CUPY and (_use_cuda() or isinstance(x, cp.ndarray))
    if use_cuda:
        x_g = x if isinstance(x, cp.ndarray) else cp.asarray(x)
        w_g = w if isinstance(w, cp.ndarray) else cp.asarray(w)
        s_g = s if isinstance(s, cp.ndarray) else cp.asarray(s)
        sim_g = sim if isinstance(sim, cp.ndarray) else cp.asarray(sim)
        sig_g = sig if isinstance(sig, cp.ndarray) else cp.asarray(sig)
        eta_g = eta if isinstance(eta, cp.ndarray) else cp.asarray(eta)

        x_f = fft2_gpu(x_g, axes=(0, 1))
        denom = ((1 + eta_g) * cp.sum(cp.abs(x_f) ** 2, axis=3, keepdims=True)) + sig_g
        c = cp.conj(x_f) / denom

        w_f = fft2_gpu(w_g, axes=(0, 1))
        w_ff = fft2_gpu(w_f, axes=(0, 1))
        s_term = s_f if s_f is not None else fft2_gpu(s_g, axes=(0, 1))
        sim_term = sim_f if sim_f is not None else fft2_gpu(sim_g, axes=(0, 1))
        rf = s_term + eta_g * sim_term - cp.sum(x_f * w_ff, axis=3, keepdims=True)
        w_f_rep = cp.broadcast_to(w_f, (w_f.shape[0], w_f.shape[1], para.ch, w_f.shape[3], w_f.shape[4]))
        gf = w_f_rep + c * rf

        g = cp.real(ifft2_gpu(gf, axes=(0, 1)))
        g = cp.sum(g, axis=2, keepdims=True) / para.ch
        if return_debug:
            dbg = {"C": c, "Rf": rf, "Gf": gf}
            return g, dbg
        return g

    x_f = _fft2(x)
    denom = ((1 + eta) * np.sum(np.abs(x_f) ** 2, axis=3, keepdims=True)) + sig
    c = np.conj(x_f) / denom

    w_f = _fft2(w)
    s_term = s_f if s_f is not None else _fft2(s)
    sim_term = sim_f if sim_f is not None else _fft2(sim)
    rf = s_term + eta * sim_term - np.sum(x_f * _fft2(w_f), axis=3, keepdims=True)
    w_f_rep = np.broadcast_to(w_f, (w_f.shape[0], w_f.shape[1], para.ch, w_f.shape[3], w_f.shape[4]))
    gf = w_f_rep + c * rf

    g = _ifft2(gf)
    g = np.sum(g, axis=2, keepdims=True) / para.ch
    if return_debug:
        dbg = {"C": c, "Rf": rf, "Gf": gf}
        return g, dbg
    return g
