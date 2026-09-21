import numpy as np
from scipy import ndimage, signal
from utils.matlab_engine_fft import engine_available, ifft2_symmetric_via_engine
from utils.fft_backend import ifft2 as ifft2_backend


def im2double(img):
    if np.issubdtype(img.dtype, np.floating):
        return img.astype(np.float64)
    if np.issubdtype(img.dtype, np.integer):
        info = np.iinfo(img.dtype)
        return img.astype(np.float64) / info.max
    return img.astype(np.float64)


def im2uint8(img):
    if np.issubdtype(img.dtype, np.floating):
        out = np.clip(img, 0.0, 1.0) * 255.0
        return np.round(out).astype(np.uint8)
    if np.issubdtype(img.dtype, np.integer):
        if img.dtype == np.uint8:
            return img
        info = np.iinfo(img.dtype)
        out = (img.astype(np.float64) / info.max) * 255.0
        return np.round(out).astype(np.uint8)
    return np.clip(img, 0, 255).astype(np.uint8)


def rgb2gray(img):
    if img.ndim == 2:
        return img
    if img.shape[2] == 1:
        return img[:, :, 0]
    r = img[:, :, 0]
    g = img[:, :, 1]
    b = img[:, :, 2]
    return 0.2989 * r + 0.5870 * g + 0.1140 * b


def gray2ind(gray, n):
    gray = np.clip(gray, 0.0, 1.0)
    idx = np.floor(gray * (n - 1)).astype(np.int64) + 1
    idx = np.clip(idx, 1, n)
    return idx


def ind2rgb(ind, cmap):
    flat = ind.reshape(-1).astype(np.int64)
    flat = np.clip(flat, 1, cmap.shape[0]) - 1
    rgb = cmap[flat]
    return rgb.reshape(ind.shape + (3,))


def jet(m=256):
    if m <= 0:
        return np.zeros((0, 3), dtype=np.float64)
    x = np.linspace(0.0, 1.0, m)
    r = np.clip(1.5 - np.abs(4.0 * x - 3.0), 0.0, 1.0)
    g = np.clip(1.5 - np.abs(4.0 * x - 2.0), 0.0, 1.0)
    b = np.clip(1.5 - np.abs(4.0 * x - 1.0), 0.0, 1.0)
    return np.stack([r, g, b], axis=1)


def _mode_map(mode):
    if mode == "replicate":
        return "nearest"
    if mode == "symmetric":
        return "reflect"
    if mode == "circular":
        return "wrap"
    return "constant"


def _normalize_kernel(kernel, ndim):
    k = np.asarray(kernel, dtype=np.float64)
    if k.ndim == 1 and ndim >= 2:
        k = k.reshape(1, -1)
    return k


def imfilter(img, kernel, mode="replicate"):
    mode_mapped = _mode_map(mode)
    k = _normalize_kernel(kernel, img.ndim)
    if img.ndim == 2:
        return ndimage.correlate(img, k, mode=mode_mapped)
    if img.ndim == 3:
        out = np.zeros_like(img, dtype=np.float64)
        for c in range(img.shape[2]):
            out[:, :, c] = ndimage.correlate(img[:, :, c], k, mode=mode_mapped)
        return out
    raise ValueError("imfilter expects 2D or 3D array")


def conv2(a, b, mode="full"):
    if a.ndim == 1 or b.ndim == 1:
        return np.convolve(a.ravel(), b.ravel(), mode=mode)
    return signal.convolve2d(a, b, mode=mode, boundary="fill", fillvalue=0)


def padarray(arr, pad_width, mode="constant", direction="both"):
    if mode == "post":
        direction = "post"
        mode = "constant"
    if isinstance(pad_width, (list, tuple)):
        pad_width = list(pad_width)
    else:
        pad_width = [pad_width]
    if len(pad_width) == 2:
        pad_width = [pad_width[0], pad_width[1]] + ([0] * (arr.ndim - 2))
    if len(pad_width) < arr.ndim:
        pad_width = pad_width + [0] * (arr.ndim - len(pad_width))

    if direction == "post":
        pad = [(0, p) for p in pad_width]
    elif direction == "pre":
        pad = [(p, 0) for p in pad_width]
    else:
        pad = [(p, p) for p in pad_width]

    if mode == "symmetric":
        return np.pad(arr, pad, mode="symmetric")
    if mode == "replicate":
        return np.pad(arr, pad, mode="edge")
    return np.pad(arr, pad, mode="constant")


def fspecial_gaussian(shape, sigma):
    if isinstance(shape, (list, tuple)):
        rows, cols = shape
    else:
        rows = cols = shape
    r = np.arange(-(rows // 2), rows // 2 + 1)
    c = np.arange(-(cols // 2), cols // 2 + 1)
    rr, cc = np.meshgrid(r, c, indexing="ij")
    h = np.exp(-(rr ** 2 + cc ** 2) / (2.0 * sigma ** 2))
    h /= np.sum(h)
    return h


def ifft2_symmetric(x):
    """
    MATLAB-compatible ifft2(...,'symmetric') over the first two axes.
    Enforces Hermitian symmetry in frequency before inverse FFT.
    """
    if engine_available():
        return ifft2_symmetric_via_engine(x)

    x = np.asarray(x)
    m, n = x.shape[0], x.shape[1]
    ii = (-np.arange(m)) % m
    jj = (-np.arange(n)) % n
    x_pair = np.conj(x[np.ix_(ii, jj)])
    x_sym = 0.5 * (x + x_pair)

    # Frequencies that map to themselves must be real-valued.
    i_self = [0]
    j_self = [0]
    if m % 2 == 0:
        i_self.append(m // 2)
    if n % 2 == 0:
        j_self.append(n // 2)
    for i in i_self:
        for j in j_self:
            x_sym[i, j, ...] = np.real(x_sym[i, j, ...])

    y = ifft2_backend(x_sym, axes=(0, 1))
    return np.real(y)
