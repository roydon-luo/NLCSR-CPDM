import os
import numpy as np
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


def Z_update(w, d, s, sim, rho, eta, para, return_debug=False, s_f=None, sim_f=None):
    use_cuda = _HAS_CUPY and (_use_cuda() or isinstance(w, cp.ndarray))
    if use_cuda:
        w_g = w if isinstance(w, cp.ndarray) else cp.asarray(w)
        d_g = d if isinstance(d, cp.ndarray) else cp.asarray(d)
        s_g = s if isinstance(s, cp.ndarray) else cp.asarray(s)
        sim_g = sim if isinstance(sim, cp.ndarray) else cp.asarray(sim)
        eta_g = eta if isinstance(eta, cp.ndarray) else cp.asarray(eta)
        rho_g = rho if isinstance(rho, cp.ndarray) else cp.asarray(rho)

        d_f = fft2_gpu(d_g, axes=(0, 1))[..., cp.newaxis]
        denom = ((1 + eta_g) * cp.sum(cp.abs(d_f) ** 2, axis=3, keepdims=True)) + rho_g
        c = cp.conj(d_f) / denom

        w_f = fft2_gpu(w_g, axes=(0, 1))
        s_term = s_f if s_f is not None else fft2_gpu(s_g, axes=(0, 1))
        sim_term = sim_f if sim_f is not None else fft2_gpu(sim_g, axes=(0, 1))
        r = s_term + eta_g * sim_term - (1 + eta_g) * cp.sum(w_f * d_f, axis=3, keepdims=True)
        zf = w_f + c * r
        z = cp.real(ifft2_gpu(zf, axes=(0, 1)))
        if return_debug:
            return z, {"C": c, "R": r, "Zf": zf}
        return z

    d_f = _fft2(d)[..., np.newaxis]
    denom = ((1 + eta) * np.sum(np.abs(d_f) ** 2, axis=3, keepdims=True)) + rho
    c = np.conj(d_f) / denom

    w_f = _fft2(w)
    s_term = s_f if s_f is not None else _fft2(s)
    sim_term = sim_f if sim_f is not None else _fft2(sim)
    r = s_term + eta * sim_term - (1 + eta) * np.sum(w_f * d_f, axis=3, keepdims=True)
    zf = w_f + c * r

    z = _ifft2(zf)
    if return_debug:
        return z, {"C": c, "R": r, "Zf": zf}
    return z
