import os

import numpy as np

_ENG = None
_ENGINE_READY = None


def _use_engine_reduction():
    return os.environ.get("USE_MATLAB_ENGINE_REDUCTION", "").strip() == "1"


def _start_engine():
    global _ENG, _ENGINE_READY
    if _ENGINE_READY is not None:
        return _ENGINE_READY
    try:
        import matlab.engine  # type: ignore

        _ENG = matlab.engine.start_matlab("-nojvm")
        _ENGINE_READY = True
    except Exception:
        _ENG = None
        _ENGINE_READY = False
    return _ENGINE_READY


def enabled():
    if not _use_engine_reduction():
        return False
    return _start_engine()


def _to_matlab_complex(arr):
    import matlab  # type: ignore

    arr = np.asarray(arr)
    if np.iscomplexobj(arr):
        return matlab.double(arr.tolist(), is_complex=True)
    return matlab.double(arr.tolist())


def sum_abs2(arr):
    if not enabled():
        return float(np.sum(np.abs(np.asarray(arr).ravel()) ** 2))
    m_arr = _to_matlab_complex(np.asarray(arr))
    _ENG.workspace["arr_py"] = m_arr
    _ENG.eval("rho_py = sum(abs(arr_py(:)).^2);", nargout=0)
    rho = float(_ENG.workspace["rho_py"])
    _ENG.eval("clear arr_py rho_py;", nargout=0)
    return rho


def sum_conj_mul(p, q):
    if not enabled():
        return float(np.sum((np.conj(np.asarray(p)) * np.asarray(q)).ravel()).real)
    m_p = _to_matlab_complex(np.asarray(p))
    m_q = _to_matlab_complex(np.asarray(q))
    _ENG.workspace["p_py"] = m_p
    _ENG.workspace["q_py"] = m_q
    _ENG.eval("pq_py = sum(reshape(conj(p_py).*q_py,[],1));", nargout=0)
    pq = _ENG.workspace["pq_py"]
    _ENG.eval("clear p_py q_py pq_py;", nargout=0)
    return float(np.real(pq))

