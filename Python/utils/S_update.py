import os
import numpy as np
from scipy import signal
from utils.pcg_ND import pcg_ND

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False


def S_update(dx, mask, para, return_debug=False):
    h = para.h
    w = para.w
    ch = para.ch
    beta_hub = para.beta_hub
    rho_hub = para.rho_hub
    delta_hub = para.delta_hub
    maxiter = para.maxiter
    tol = para.tol
    s_mask = para.Smask
    gamma = para.gamma

    use_cuda = _HAS_CUPY and isinstance(dx, cp.ndarray) and os.environ.get("USE_CUDA", "0").strip() == "1"
    xp = cp if use_cuda else np

    dx = xp.asarray(dx, dtype=xp.float64)
    mask = xp.asarray(mask, dtype=xp.float64)
    s_mask = xp.asarray(s_mask, dtype=xp.float64)
    beta_hub = xp.asarray(beta_hub, dtype=xp.float64)
    rho_hub = xp.asarray(rho_hub, dtype=xp.float64)
    delta_hub = xp.asarray(delta_hub, dtype=xp.float64)

    def norm2(x):
        return float(xp.linalg.norm(x.ravel()))

    def rho_hub_update(x):
        rh = xp.reshape(rho_hub, (1, 1, ch))
        return x * rh

    def huber(x, delta, rho, lamb):
        eps = np.finfo(np.float64).eps
        t = (lamb + rho) / (rho + eps) * delta
        out = xp.zeros_like(x)
        mask1 = xp.abs(x) <= t
        mask2 = x > t
        mask3 = x < -t
        out[mask1] = (rho / (lamb + rho + eps)) * x[mask1]
        out[mask2] = x[mask2] - (delta / (rho + eps)) * lamb
        out[mask3] = x[mask3] + (delta / (rho + eps)) * lamb
        return out

    dx_k_cpu = np.array([[0, -1, 1]], dtype=np.float64)
    dy_k_cpu = dx_k_cpu.T
    dxx_cpu = signal.convolve2d(dx_k_cpu, dx_k_cpu, mode="full", boundary="fill", fillvalue=0)
    dyy_cpu = signal.convolve2d(dy_k_cpu, dy_k_cpu, mode="full", boundary="fill", fillvalue=0)
    dxy_cpu = signal.convolve2d(dx_k_cpu, dy_k_cpu, mode="full", boundary="fill", fillvalue=0)

    dx_k = xp.asarray(dx_k_cpu, dtype=xp.float64)
    dy_k = xp.asarray(dy_k_cpu, dtype=xp.float64)
    dxx = xp.asarray(dxx_cpu, dtype=xp.float64)
    dyy = xp.asarray(dyy_cpu, dtype=xp.float64)
    dxy = xp.asarray(dxy_cpu, dtype=xp.float64)

    def _correlate_circular_2d(img2, ker2):
        kh, kw = int(ker2.shape[0]), int(ker2.shape[1])
        ch0, cw0 = kh // 2, kw // 2
        out2 = xp.zeros_like(img2, dtype=xp.float64)
        for u in range(kh):
            du = u - ch0
            for v in range(kw):
                c0 = ker2[u, v]
                if float(c0) == 0.0:
                    continue
                dv = v - cw0
                out2 = out2 + c0 * xp.roll(img2, shift=(-du, -dv), axis=(0, 1))
        return out2

    def imfilter_wrap(img, kernel):
        if img.ndim == 2:
            return _correlate_circular_2d(img, kernel)
        out = xp.zeros_like(img, dtype=xp.float64)
        for cc in range(img.shape[2]):
            out[:, :, cc] = _correlate_circular_2d(img[:, :, cc], kernel)
        return out

    def Dx(x):
        return imfilter_wrap(x, dx_k)

    def Dy(x):
        return imfilter_wrap(x, dy_k)

    def Dxx(x):
        return imfilter_wrap(x, dxx)

    def Dyy(x):
        return imfilter_wrap(x, dyy)

    def Dxy(x):
        return imfilter_wrap(x, dxy)

    def Diff(v):
        return xp.stack([Dx(v), Dy(v), Dxx(v), Dxy(v), Dyy(v)], axis=3)

    def DxT(x):
        return imfilter_wrap(x, xp.rot90(dx_k, 2))

    def DyT(x):
        return imfilter_wrap(x, xp.rot90(dy_k, 2))

    def DxxT(x):
        return imfilter_wrap(x, xp.rot90(dxx, 2))

    def DyyT(x):
        return imfilter_wrap(x, xp.rot90(dyy, 2))

    def DxyT(x):
        return imfilter_wrap(x, xp.rot90(dxy, 2))

    def DiffT(w_in):
        return DxT(w_in[:, :, :, 0]) + DyT(w_in[:, :, :, 1]) + DxxT(w_in[:, :, :, 2]) + DxyT(w_in[:, :, :, 3]) + DyyT(
            w_in[:, :, :, 4]
        )

    s = xp.asarray(
        [
            [0, 0, 1, 0, 0],
            [0, 1, -7, 1, 0],
            [1, -7, 20, -7, 1],
            [0, 1, -7, 1, 0],
            [0, 0, 1, 0, 0],
        ],
        dtype=xp.float64,
    )

    def DTD(v):
        return imfilter_wrap(v, s)

    gamma_fix = float(gamma)
    b1 = dx + gamma_fix * s_mask

    def Fi(x):
        return (gamma_fix * mask + 1.0) * x

    def A(v):
        return Fi(v) + rho_hub_update(DTD(v))

    z_hub = xp.zeros((h, w, ch, 5), dtype=xp.float64)
    u_hub = xp.zeros((h, w, ch, 5), dtype=xp.float64)
    obj_total = np.zeros((maxiter + 2, 3), dtype=np.float64)
    obj_total[0, :] = np.inf
    dbg = {}
    if return_debug:
        dbg["b_hist"] = []
        dbg["x_hist"] = []
        dbg["gamma_hist"] = []
        dbg["pcg_hist"] = []

    itr = 0
    x = xp.zeros_like(dx)
    while itr <= maxiter:
        itr += 1
        b2 = rho_hub_update(DiffT(z_hub - u_hub))
        b = b1 + b2
        if return_debug:
            x, dbg_pcg = pcg_ND(A, b, tol, return_debug=True)
        else:
            x = pcg_ND(A, b, tol)
        d = Diff(x)
        d_u_hub = d + u_hub
        for ch_hub in range(ch):
            z_hub[:, :, ch_hub, :] = huber(d_u_hub[:, :, ch_hub, :], delta_hub[ch_hub], rho_hub[ch_hub], beta_hub[ch_hub])
        u_hub = d_u_hub - z_hub
        if return_debug and itr == 1:
            dbg["itr1_b1"] = b1.copy()
            dbg["itr1_b2"] = b2.copy()
            dbg["itr1_b"] = b.copy()
            dbg["itr1_x"] = x.copy()
            dbg["itr1_d"] = d.copy()
            dbg["itr1_Z_hub"] = z_hub.copy()
            dbg["itr1_U_hub"] = u_hub.copy()
            dbg["itr1_pcg"] = dbg_pcg
        if return_debug:
            dbg["b_hist"].append(b.copy())
            dbg["x_hist"].append(x.copy())
            dbg["gamma_hist"].append(float(gamma))
            dbg["pcg_hist"].append(dbg_pcg)

        obj_data = 0.5 * (norm2(x - b1) ** 2 + gamma_fix * norm2(s_mask - mask * x) ** 2)
        obj_hub = 0.0
        for ch_hub in range(ch):
            obj_hub += beta_hub[ch_hub] * objective_hub(d[:, :, ch_hub, :], delta_hub[ch_hub], xp)
        obj_total[itr, 0] = obj_data + float(obj_hub)
        obj_total[itr, 1] = obj_data
        obj_total[itr, 2] = float(obj_hub)
        rel_err = np.abs(obj_total[itr, 0] - obj_total[itr - 1, 0])
        if rel_err < tol:
            break
        gamma *= 0.9

    if return_debug:
        dbg["final_x"] = x.copy()
        dbg["final_itr"] = itr
        if len(dbg["b_hist"]):
            dbg["b_hist"] = xp.stack(dbg["b_hist"], axis=3)
            dbg["x_hist"] = xp.stack(dbg["x_hist"], axis=3)
        else:
            dbg["b_hist"] = xp.zeros((h, w, ch, 0), dtype=xp.float64)
            dbg["x_hist"] = xp.zeros((h, w, ch, 0), dtype=xp.float64)
        dbg["gamma_hist"] = np.asarray(dbg["gamma_hist"], dtype=np.float64)
        return x, dbg
    return x


def objective_hub(z, delta, xp):
    return xp.sum(huber_loss(z.ravel(), delta, xp))


def huber_loss(z, k, xp):
    g_hub = xp.zeros_like(z, dtype=xp.float64)
    if float(k) != 0.0:
        mask1 = xp.abs(z) <= k
        mask2 = xp.abs(z) > k
        g_hub[mask1] = 0.5 * (z[mask1] ** 2)
        g_hub[mask2] = k * (xp.abs(z[mask2]) - 0.5 * k)
    else:
        g_hub = xp.abs(z)
    return g_hub
