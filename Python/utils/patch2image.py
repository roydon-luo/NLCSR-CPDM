import numpy as np

try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False


def patch2image(s_recblk, para):
    """
    Vectorized patch2image reconstruction.
    Optimized for GPU (CuPy) to avoid slow Python loops.
    """
    use_cuda = _HAS_CUPY and isinstance(s_recblk, cp.ndarray)
    xp = cp if use_cuda else np

    f = para.f
    ch = para.ch
    h = para.h
    w = para.w

    blk = xp.zeros((h, w, ch), dtype=xp.float64)
    cnt = xp.zeros((h, w, ch), dtype=xp.float64)
    
    # Ensure inputs are correct type/device
    s_recblk = xp.asarray(s_recblk, dtype=xp.float64)

    # Convert para lists to arrays
    r = xp.asarray(para.r)
    c = xp.asarray(para.c)
    
    n_row = len(r)
    n_col = len(c)
    n_patches_expected = n_row * n_col
    
    # Validation (optional but recommended)
    if s_recblk.shape[3] != n_patches_expected:
        # If sizes mismatch, fallback or allow broadcasting if designed so.
        # Assuming strictly matching dimensions for this optimization.
        pass

    # --- Vectorized Coordinate Generation ---
    # Original loop order: for j in c: for i in r:
    # Outer loop: cols (c), Inner loop: rows (r)
    # We construct a grid for all top-left corners.
    
    # J repeats c for every full traversal of r (c0, c0..., c1, c1...)
    J_grid = xp.repeat(c, n_row) 
    # I tiles r for every c (r0, r1..., r0, r1...)
    I_grid = xp.tile(r, n_col)
    
    # Generate pixel offsets within a patch (f x f)
    # dy, dx shape: (f, f)
    dy, dx = xp.indices((f, f))
    
    # Broadcast to (N_patches, f, f)
    # I_grid: (N,) -> (N, 1, 1)
    Y = I_grid[:, None, None] + dy[None, :, :]
    X = J_grid[:, None, None] + dx[None, :, :]
    
    # Flatten spatial coordinates to 1D
    Y_flat = Y.ravel()
    X_flat = X.ravel()
    
    # Prepare patch values
    # s_recblk is (f, f, ch, N). We need to align it with Y, X (which are N-major)
    # Move N to front: (N, f, f, ch)
    patches = xp.moveaxis(s_recblk, 3, 0)
    # Flatten spatial dimensions: (N*f*f, ch)
    patches_flat = patches.reshape(-1, ch)
    
    # --- Atomic Accumulation ---
    # xp.add.at works for both numpy and cupy
    # We use (Y_flat, X_flat) as indices. 
    # The last dimension (ch) is broadcasted automatically.
    xp.add.at(blk, (Y_flat, X_flat), patches_flat)
    xp.add.at(cnt, (Y_flat, X_flat), 1.0)

    s_rec = blk / xp.maximum(cnt, xp.finfo(xp.float64).eps)
    return s_rec