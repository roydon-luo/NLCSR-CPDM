import numpy as np
from utils.matlab_compat import rgb2gray, gray2ind, ind2rgb, jet, im2uint8


def colorjetmap(i):
    i_gray = rgb2gray(i)
    cmap = jet(256)
    ind = gray2ind(i_gray, 256)
    i_jet = ind2rgb(ind, cmap)
    return im2uint8(i_jet)
