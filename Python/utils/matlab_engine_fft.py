import os
import subprocess
import uuid
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

_ENG = None
_ENGINE_READY = None


def _use_engine():
    return os.environ.get("USE_MATLAB_ENGINE_IFFT", "").strip() == "1"


def _conda_env_name():
    return os.environ.get("MATLAB_ENGINE_CONDA_ENV", "").strip()


def _start_engine_inprocess():
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


def _worker_script_path():
    return Path(__file__).resolve().parent.parent / "matlab_bridge" / "engine_ifft2_worker.py"


def _run_worker(in_path, out_path):
    env_name = _conda_env_name()
    if not env_name:
        return False
    cmd = [
        "D:\\ProgramData\\anaconda3\\Scripts\\conda.exe",
        "run",
        "-n",
        env_name,
        "python",
        str(_worker_script_path()),
        str(in_path),
        str(out_path),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return r.returncode == 0
    except Exception:
        return False


def engine_available():
    if not _use_engine():
        return False
    if _start_engine_inprocess():
        return True
    return _conda_env_name() != "" and _worker_script_path().exists()


def _bridge_tmp_dir():
    d = Path(__file__).resolve().parent.parent / ".engine_tmp"
    d.mkdir(exist_ok=True)
    return d


def ifft2_symmetric_via_engine(x):
    x = np.asarray(x)
    td = _bridge_tmp_dir()
    uid = uuid.uuid4().hex
    in_path = td / f"in_{uid}.mat"
    out_path = td / f"out_{uid}.mat"

    try:
        savemat(in_path, {"X": x}, do_compression=False)

        if _start_engine_inprocess():
            in_posix = in_path.as_posix().replace("'", "''")
            out_posix = out_path.as_posix().replace("'", "''")
            cmd = (
                f"s=load('{in_posix}','X'); "
                "Y=ifft2(s.X,'symmetric'); "
                f"save('{out_posix}','Y','-v7');"
            )
            _ENG.eval(cmd, nargout=0)
        else:
            ok = _run_worker(in_path, out_path)
            if not ok:
                raise RuntimeError("MATLAB Engine bridge unavailable")

        y = loadmat(out_path)["Y"]
        return np.asarray(y, dtype=np.float64)
    finally:
        try:
            if in_path.exists():
                in_path.unlink()
        except Exception:
            pass
        try:
            if out_path.exists():
                out_path.unlink()
        except Exception:
            pass
