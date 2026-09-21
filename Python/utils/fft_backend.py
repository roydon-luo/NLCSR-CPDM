import os
import numpy as np

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False

# DFT-only backend. No FFT path is used.
_DFT_CACHE_NP = {}
_DFT_CACHE_CP = {}


def _use_cuda():
    return os.environ.get('USE_CUDA', '0').strip() == '1' and _HAS_CUPY


def _axes(axes, ndim):
    if len(axes) != 2:
        raise ValueError('Only 2D transform is supported')
    a0 = axes[0] % ndim
    a1 = axes[1] % ndim
    if a0 == a1:
        raise ValueError('Invalid axes')
    return a0, a1


def _dft_np(n):
    n = int(n)
    if n in _DFT_CACHE_NP:
        return _DFT_CACHE_NP[n]
    k = np.arange(n, dtype=np.float64).reshape(-1, 1)
    m = np.arange(n, dtype=np.float64).reshape(1, -1)
    w = np.exp(-2j * np.pi * (k @ m) / n)
    iw = np.conj(w) / n
    _DFT_CACHE_NP[n] = (w.astype(np.complex128), iw.astype(np.complex128))
    return _DFT_CACHE_NP[n]


def _dft_cp(n):
    n = int(n)
    if n in _DFT_CACHE_CP:
        return _DFT_CACHE_CP[n]
    k = cp.arange(n, dtype=cp.float64).reshape(-1, 1)
    m = cp.arange(n, dtype=cp.float64).reshape(1, -1)
    w = cp.exp(-2j * cp.pi * (k @ m) / n)
    iw = cp.conj(w) / n
    _DFT_CACHE_CP[n] = (w.astype(cp.complex128), iw.astype(cp.complex128))
    return _DFT_CACHE_CP[n]


def _dft2_np(x, axes=(0, 1), inverse=False):
    a0, a1 = _axes(axes, x.ndim)
    z = np.moveaxis(x, (a0, a1), (0, 1)).astype(np.complex128, copy=False)
    n0, n1 = z.shape[0], z.shape[1]
    w0, iw0 = _dft_np(n0)
    w1, iw1 = _dft_np(n1)
    a_mat0 = iw0 if inverse else w0
    a_mat1 = iw1 if inverse else w1
    t = np.tensordot(a_mat0, z, axes=(1, 0))
    y = np.tensordot(t, a_mat1.T, axes=(1, 0))
    y = np.moveaxis(y, -1, 1)
    return np.moveaxis(y, (0, 1), (a0, a1))


def _dft2_cp(x, axes=(0, 1), inverse=False):
    a0, a1 = _axes(axes, x.ndim)
    z = cp.asarray(cp.moveaxis(x, (a0, a1), (0, 1)), dtype=cp.complex128)
    n0, n1 = z.shape[0], z.shape[1]
    w0, iw0 = _dft_cp(n0)
    w1, iw1 = _dft_cp(n1)
    a_mat0 = iw0 if inverse else w0
    a_mat1 = iw1 if inverse else w1
    t = cp.tensordot(a_mat0, z, axes=(1, 0))
    y = cp.tensordot(t, a_mat1.T, axes=(1, 0))
    y = cp.moveaxis(y, -1, 1)
    y = cp.moveaxis(y, (0, 1), (a0, a1))
    return y


def _is_cupy_array(x):
    return _HAS_CUPY and isinstance(x, cp.ndarray)


def fft2_gpu(x, axes=(0, 1)):
    if not _HAS_CUPY:
        raise RuntimeError('CuPy is unavailable')
    return _dft2_cp(cp.asarray(x), axes=axes, inverse=False)


def ifft2_gpu(x, axes=(0, 1)):
    if not _HAS_CUPY:
        raise RuntimeError('CuPy is unavailable')
    return _dft2_cp(cp.asarray(x), axes=axes, inverse=True)


def fft2(x, axes=(0, 1)):
    if _is_cupy_array(x):
        return _dft2_cp(x, axes=axes, inverse=False)
    x = np.asarray(x)
    return _dft2_np(x, axes=axes, inverse=False)


def ifft2(x, axes=(0, 1)):
    if _is_cupy_array(x):
        return _dft2_cp(x, axes=axes, inverse=True)
    x = np.asarray(x)
    return _dft2_np(x, axes=axes, inverse=True)
