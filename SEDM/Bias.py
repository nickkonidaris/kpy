import numpy as np
import pyfits as pf
# Bias Subtraction


def full_frame(dat):

    bias = np.median(dat[2045:,:], axis=0)

    return np.tile(bias, (2048, 1))

def remove(fits_obj):
    '''Return the bias-subtracted version of the fits object'''

    dat = fits_obj[0].data

    if dat.shape == (2048, 2048):
        bias_img = full_frame(dat)

    return dat - bias_img
