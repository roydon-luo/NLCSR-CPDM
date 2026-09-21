import os
import numpy as np
import imageio.v2 as imageio
from utils.matlab_compat import im2double


def load_cpfa(filename):
    img = imageio.imread(filename)
    img = im2double(img)
    # Default: keep original resolution for arbitrary-size input.
    # Optional crop can be enabled by setting both env vars:
    # CPFA_MAX_H=<int>, CPFA_MAX_W=<int>
    max_h = int(os.environ.get("CPFA_MAX_H", "0"))
    max_w = int(os.environ.get("CPFA_MAX_W", "0"))
    if max_h <= 0 or max_w <= 0:
        return img

    height, width = img.shape[0], img.shape[1]
    if height <= max_h and width <= max_w:
        return img

    # Deterministic centered crop when optional cap is enabled.
    start_row = max((height - max_h) // 2, 0)
    start_col = max((width - max_w) // 2, 0)
    return img[start_row:start_row + max_h, start_col:start_col + max_w, ...]
