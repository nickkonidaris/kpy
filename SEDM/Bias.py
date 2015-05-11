import numpy as np
import glob
import pyfits as pf
import scipy.ndimage.filters as FI
# Bias Subtraction


def full_frame(dat):

    bias = np.median(dat[:,2045:], axis=1)
    bias = bias.astype(np.float)
    smooth = FI.median_filter(bias, size=50)

    return np.tile(smooth, (2048,1))

def remove(fits_obj):
    '''Return the bias-subtracted version of the fits object'''

    dat = fits_obj[0].data
    try: GAIN = fits_obj[0].header['GAIN']
    except: GAIN = 1.8 # Guess the gain

    if dat.shape == (2048, 2048):
        bias_img = full_frame(dat)

    return (dat - bias_img.T) * GAIN

def add_prefix(fname):
    '''/path/to/file --> /path/to/b_file'''

    sp = fname.split("/")
    sp[-1] = 'b_' + sp[-1]

    return "/".join(sp)

if __name__ == '__main__':
    import sys

    files = sys.argv[1:]
    
    for file in files:
        try:
            if file[-5:] != '.fits': continue
        except:
            continue
        print file
        FF = pf.open(file)
        adcspeed = FF[0].header['ADCSPEED']

        bfname = "bias%1.1f.fits" % adcspeed
        bias = pf.open(bfname)

        FF[0].data -= bias[0].data
        FF[0].data = remove(FF)

        outname = add_prefix(file)
        FF[0].header['BIASSUB'] = ('Subtracted', 'Ovrscn + bias handled by Bias.py')
        FF[0].header['BIASSUB2'] = (bfname , 'Bias file used')
        try: 
            GAIN = FF[0].header['GAIN']
            FF[0].header['GAIN'] = (1.0, 'GAIN Adjusted (was %s)' % GAIN)
        except: 
            GAIN=1.8 # Guess the gain 
            FF[0].header['GAIN'] = (1.0, 'GAIN Adjusted (was guessed %s)' % GAIN)
        FF[0].header['BUNIT'] = ('electron')
        FF.writeto(outname)

