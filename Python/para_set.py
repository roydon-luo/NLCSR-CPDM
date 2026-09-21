import os
from types import SimpleNamespace
import numpy as np

# Default to explicit DFT backend unless user overrides it externally.
os.environ.setdefault("USE_EXPLICIT_DFT", "1")


def para_set(s):
    h, w, ch = s.shape
    para = SimpleNamespace()
    para.h = h
    para.w = w
    para.ch = ch
    para.m = 5
    para.K = 16
    para.f = 6
    para.ns = 10
    para.step = 3
    para.nSig = 0
    para.C = 0.28
    para.C_min = 1e-4
    para.min_dis = 0.1
    para.Eps = 1e-6
    para.relaxParam_d = 1.8
    para.sigma = 1
    para.AutoSigma = 1
    para.SigmaUpdateCycle = 1
    para.cdl_iters = 1
    para.relaxParam_x = 1.8
    para.rho = 1
    para.AutoRho = 1
    para.RhoUpdateCycle = 1
    para.csc_iters = 1
    para.lamb = 1e-4
    para.relaxParam_s = 1.8
    para.gamma = 0.1
    para.tol = 1e-5
    para.rho_hub = np.repeat(0.05, 12)
    para.delta_hub = np.repeat(0.001, 12)
    para.beta_hub = np.repeat(0.95, 12)
    para.main_iternum = 1
    para.mu = 5.0
    para.itar = 1.2
    para.maxiter = 20
    para.eAbs = 1e-4
    para.eRel = 1e-4
    # Optimized for 24GB VRAM: Process all patches in one batch if possible
    para.patch_num = 10000 
    para.gpu_num = 0
    return para