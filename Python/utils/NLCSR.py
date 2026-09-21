import numpy as np
import os
from utils.get_simblk import get_simblk
from utils.get_simpara import get_simpara
from utils.Z_update import Z_update
from utils.G_update import G_update
from utils.matlab_compat import padarray
from utils.fft_backend import fft2 as fft2_backend, ifft2 as ifft2_backend

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False


def _fft2(x):
    return fft2_backend(x, axes=(0, 1))


def _ifft2(x):
    return ifft2_backend(x, axes=(0, 1))


def _use_cuda():
    return _HAS_CUPY and (os.environ.get("USE_CUDA", "0").strip() == "1")


def NLCSR(d0, s_pat, para, return_debug=False, return_numpy=True):
    maxiter = para.maxiter
    e_abs = para.eAbs
    e_rel = para.eRel
    ch = para.ch
    m = para.m
    k = para.K
    mu = para.mu
    itar = para.itar
    alphad = para.relaxParam_d
    sig_update_cycle = para.SigmaUpdateCycle
    cdl_iters = para.cdl_iters
    rho = para.rho
    f = para.f
    c = para.C
    alphax = para.relaxParam_x
    lamb = para.lamb
    rho_update_cycle = para.RhoUpdateCycle
    csc_iters = para.csc_iters
    sig = para.sigma

    l = s_pat.shape[3]
    sim_pat, pat_arr = get_simblk(s_pat, para)
    b = get_simpara(s_pat, sim_pat, pat_arr, para)

    use_cuda = _use_cuda() or (_HAS_CUPY and (isinstance(s_pat, cp.ndarray) or isinstance(d0, cp.ndarray)))
    xp = cp if use_cuda else np

    b = xp.reshape(b, (f, f, ch, 1, l), order="F")
    d = xp.reshape(d0, (m, m, 1, k), order="F")
    d = padarray(d, [f - m, f - m], mode="post")
    s_pat = xp.reshape(s_pat, (f, f, ch, 1, l), order="F")
    sim_pat = xp.reshape(sim_pat, (f, f, ch, 1, l), order="F")
    if use_cuda:
        d = cp.asarray(d)
        s_pat = cp.asarray(s_pat)
        sim_pat = cp.asarray(sim_pat)
        b = cp.asarray(b)
    s_f = _fft2(s_pat)
    sim_f = _fft2(sim_pat)
    x = xp.zeros((f, f, ch, k, l), dtype=np.float64)
    u = xp.zeros_like(x)
    v = xp.zeros((f, f, 1, k, l), dtype=np.float64)
    nx = x.size
    nd = d.size

    itr_dx = 0
    eprix = 0
    eduax = 0
    rx = np.inf
    sx = np.inf
    eprid = 0
    eduad = 0
    rd = np.inf
    sd = np.inf
    dbg = {}

    if return_debug:
        dbg["sim_pat"] = sim_pat.copy()
        dbg["pat_arr"] = pat_arr.copy()
        dbg["B"] = b.copy()

    while itr_dx <= maxiter and (rx > eprix or sx > eduax or rd > eprid or sd > eduad):
        itr_dx += 1
        eta = c / b
        if c >= para.C_min:
            c = c * 0.92

        for tx in range(csc_iters):
            xprv = x
            w_cur = x - u
            if return_debug and itr_dx == 1 and tx == 0:
                dbg["itr1_tx1_W"] = cp.asnumpy(w_cur.copy()) if use_cuda else w_cur.copy()
            if return_debug and itr_dx == 1 and tx == 0:
                z, dbg_z = Z_update(w_cur, d, s_pat, sim_pat, rho, eta, para, return_debug=True, s_f=s_f, sim_f=sim_f)
            else:
                z = Z_update(w_cur, d, s_pat, sim_pat, rho, eta, para, s_f=s_f, sim_f=sim_f)
            if return_debug and itr_dx == 1 and tx == 0:
                dbg["itr1_tx1_Z"] = cp.asnumpy(z.copy()) if use_cuda else z.copy()
                dbg["itr1_z_C"] = cp.asnumpy(dbg_z["C"].copy()) if use_cuda else dbg_z["C"].copy()
                dbg["itr1_z_R"] = cp.asnumpy(dbg_z["R"].copy()) if use_cuda else dbg_z["R"].copy()
                dbg["itr1_z_Zf"] = cp.asnumpy(dbg_z["Zf"].copy()) if use_cuda else dbg_z["Zf"].copy()
            zr = alphax * z + (1 - alphax) * x
            x = sfthrsh(zr + u, lamb / rho)
            u = zr - x + u
            if return_debug and itr_dx == 1 and tx == 0:
                dbg["itr1_tx1_X"] = cp.asnumpy(x.copy()) if use_cuda else x.copy()
                dbg["itr1_tx1_U"] = cp.asnumpy(u.copy()) if use_cuda else u.copy()

        n_x = float(xp.linalg.norm(x.ravel()))
        n_z = float(xp.linalg.norm(z.ravel()))
        n_u = float(xp.linalg.norm(u.ravel()))
        rx = float(xp.linalg.norm((x - z).ravel()))
        sx = float(rho * xp.linalg.norm((xprv - x).ravel()))
        eprix = np.sqrt(nx) * e_abs + max(n_x, n_z) * e_rel
        eduax = np.sqrt(nx) * e_abs + rho * n_u * e_rel
        if para.AutoRho and itr_dx % rho_update_cycle == 0:
            rho, u = par_update(rho, rx, sx, mu, itar, u)
        if return_debug and itr_dx == 1:
            dbg["itr1_eta"] = cp.asnumpy(eta.copy()) if use_cuda else eta.copy()
            dbg["itr1_rho"] = rho
            dbg["itr1_Z"] = cp.asnumpy(z.copy()) if use_cuda else z.copy()
            dbg["itr1_X"] = cp.asnumpy(x.copy()) if use_cuda else x.copy()
            dbg["itr1_U"] = cp.asnumpy(u.copy()) if use_cuda else u.copy()

        for _ in range(cdl_iters):
            dprv = d
            d_exp = d[:, :, :, :, None]
            if return_debug and itr_dx == 1:
                g, dbg_g = G_update(x, d_exp - v, s_pat, sim_pat, sig, eta, para, return_debug=True, s_f=s_f, sim_f=sim_f)
            else:
                g = G_update(x, d_exp - v, s_pat, sim_pat, sig, eta, para, s_f=s_f, sim_f=sim_f)
            gr = alphad * g + (1 - alphad) * d_exp
            d = d_proj(xp.sum(gr + v, axis=4) / l)
            v = gr - d[:, :, :, :, None] + v

        n_g = float(xp.linalg.norm(g.ravel()))
        n_d = float(xp.linalg.norm(d.ravel()) * np.sqrt(l))
        n_v = float(xp.linalg.norm(v.ravel()))
        rd = float(xp.linalg.norm((g - d[:, :, :, :, None]).ravel()))
        sd = float(sig * xp.linalg.norm((dprv - d).ravel()))
        eprid = np.sqrt(nd) * e_abs + max(n_d, n_g) * e_rel
        eduad = np.sqrt(nd) * e_abs + sig * n_v * e_rel
        if para.AutoSigma and itr_dx % sig_update_cycle == 0:
            sig, v = par_update(sig, rd, sd, mu, itar, v)
        if return_debug and itr_dx == 1:
            dbg["itr1_sig"] = sig
            dbg["itr1_G"] = cp.asnumpy(g.copy()) if use_cuda else g.copy()
            dbg["itr1_D"] = cp.asnumpy(d.copy()) if use_cuda else d.copy()
            dbg["itr1_V"] = cp.asnumpy(v.copy()) if use_cuda else v.copy()
            dbg["itr1_g_C"] = cp.asnumpy(dbg_g["C"].copy()) if use_cuda else dbg_g["C"].copy()
            dbg["itr1_g_Rf"] = cp.asnumpy(dbg_g["Rf"].copy()) if use_cuda else dbg_g["Rf"].copy()
            dbg["itr1_g_Gf"] = cp.asnumpy(dbg_g["Gf"].copy()) if use_cuda else dbg_g["Gf"].copy()

    dx = _ifft2(_fft2(d)[:, :, :, :, None] * _fft2(x))
    dx = xp.real(dx)
    dx = xp.sum(dx, axis=3)
    d = xp.reshape(d[0:m, 0:m, :], (m, m, k), order="F")
    if use_cuda and return_numpy:
        dx = cp.asnumpy(dx)
        d = cp.asnumpy(d)
    if return_debug:
        dbg["final_rho"] = rho
        dbg["final_sig"] = sig
        dbg["final_C"] = c
        dbg["final_itr_dx"] = itr_dx
        return dx, d, dbg
    return dx, d


def sfthrsh(x, kappa):
    if _HAS_CUPY and isinstance(x, cp.ndarray):
        return cp.sign(x) * cp.maximum(0, cp.abs(x) - kappa)
    return np.sign(x) * np.maximum(0, np.abs(x) - kappa)


def d_proj(d):
    if _HAS_CUPY and isinstance(d, cp.ndarray):
        denom = cp.sqrt(cp.sum(d ** 2, axis=(0, 1), keepdims=True))
        denom = cp.maximum(denom, cp.finfo(cp.float64).eps)
        return d / denom
    denom = np.sqrt(np.sum(d ** 2, axis=(0, 1), keepdims=True))
    denom = np.maximum(denom, np.finfo(np.float64).eps)
    return d / denom


def par_update(par, r, s, mu, itar, y):
    a = 1
    if r > mu * s:
        a = itar
    if s > mu * r:
        a = 1 / itar
    par_ = a * par
    if par_ > 1e-4:
        par = par_
        y = y / a
    return par, y
