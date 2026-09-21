import numpy as np

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False


def _xp(arr):
    if _HAS_CUPY and isinstance(arr, cp.ndarray):
        return cp
    return np


def _sum_abs2(x):
    xp = _xp(x)
    return xp.sum(xp.abs(x) ** 2)


def _sum_conj_mul(a, b):
    xp = _xp(a)
    return xp.sum(xp.conj(a) * b)


def pcg_ND(a_fun, b, tol, maxit=100, x=None, return_debug=False):
    dbg = {
        "alpha": [],
        "beta": [],
        "rho": [],
        "pq": [],
        "normr": [],
        "stag": [],
        "moresteps": [],
        "exit_flag": "maxit_or_break",
        "iters": 0,
    }
    xp = _xp(b)
    if x is None:
        x = xp.zeros_like(b, dtype=xp.float64)
    else:
        x = xp.asarray(x, dtype=xp.float64)

    n2b = float(xp.linalg.norm(b.ravel()))
    dbg["n2b"] = float(n2b)
    if n2b == 0:
        out = xp.zeros_like(b, dtype=xp.float64)
        dbg["exit_flag"] = "zero_rhs"
        if return_debug:
            return out, dbg
        return out

    tolb = tol * n2b
    dbg["tolb"] = float(tolb)
    r = b - a_fun(x)
    normr = float(xp.linalg.norm(r.ravel()))
    if normr <= tolb:
        dbg["exit_flag"] = "initial_converged"
        if return_debug:
            return x, dbg
        return x

    rho = 1.0
    stag = 0
    moresteps = 0
    maxmsteps = min(int(np.floor(b.size / 50)), 5, b.size - maxit)
    maxstagsteps = 3

    for ii in range(maxit):
        rho1 = rho
        rho = _sum_abs2(r)
        dbg["rho"].append(float(np.real(rho)))
        if rho == 0 or np.isinf(rho):
            dbg["exit_flag"] = "rho_invalid"
            break
        if ii == 0:
            p = r.copy()
            beta = 0.0
            dbg["p_step1"] = p.copy()
        else:
            beta = rho / rho1
            if beta == 0 or np.isinf(beta):
                dbg["exit_flag"] = "beta_invalid"
                break
            p = r + beta * p
        dbg["beta"].append(float(np.real(beta)))
        q = a_fun(p)
        if ii == 0:
            dbg["q_step1"] = q.copy()
        pq = _sum_conj_mul(p, q)
        dbg["pq"].append(float(np.real(pq)))
        if pq <= 0 or np.isinf(pq):
            dbg["exit_flag"] = "pq_invalid"
            break
        alpha = rho / pq
        dbg["alpha"].append(float(np.real(alpha)))
        if np.isinf(alpha):
            dbg["exit_flag"] = "alpha_invalid"
            break

        if float(xp.linalg.norm(p.ravel())) * abs(alpha) < np.finfo(np.float64).eps * float(xp.linalg.norm(x.ravel())):
            stag += 1
        else:
            stag = 0

        x = x + alpha * p
        r = r - alpha * q
        normr = float(xp.linalg.norm(r.ravel()))
        dbg["normr"].append(float(normr))
        dbg["stag"].append(int(stag))
        dbg["moresteps"].append(int(moresteps))

        if normr <= tolb or stag >= maxstagsteps or moresteps:
            r = b - a_fun(x)
            if float(xp.linalg.norm(r.ravel())) <= tolb:
                dbg["exit_flag"] = "converged"
                break
            if stag >= maxstagsteps and moresteps == 0:
                stag = 0
            moresteps += 1
            if moresteps >= maxmsteps:
                dbg["exit_flag"] = "max_moresteps"
                break
        if stag >= maxstagsteps:
            dbg["exit_flag"] = "stagnation"
            break

    dbg["iters"] = len(dbg["rho"])
    if return_debug:
        return x, dbg
    return x
