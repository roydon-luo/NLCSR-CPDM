import numpy as np


def boxfilter(im_src, h, v):
    hei, wid = im_src.shape
    im_dst = np.zeros_like(im_src, dtype=np.float64)

    if v != 0:
        im_cum = np.cumsum(im_src, axis=0)
        im_dst[0:v + 1, :] = im_cum[v:2 * v + 1, :]
        im_dst[v + 1:hei - v, :] = im_cum[2 * v + 1:hei, :] - im_cum[0:hei - 2 * v - 1, :]
        im_dst[hei - v:hei, :] = np.tile(im_cum[hei - 1:hei, :], (v, 1)) - im_cum[hei - 2 * v - 1:hei - v - 1, :]

    if h != 0:
        if v != 0:
            im_cum = np.cumsum(im_dst, axis=1)
        else:
            im_cum = np.cumsum(im_src, axis=1)
        im_dst[:, 0:h + 1] = im_cum[:, h:2 * h + 1]
        im_dst[:, h + 1:wid - h] = im_cum[:, 2 * h + 1:wid] - im_cum[:, 0:wid - 2 * h - 1]
        im_dst[:, wid - h:wid] = np.tile(im_cum[:, wid - 1:wid], (1, h)) - im_cum[:, wid - 2 * h - 1:wid - h - 1]

    return im_dst
