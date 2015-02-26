
import argparse
import pdb
import numpy as np
import pylab as pl
import pyfits as pf
import sys


import NPK.Fit as FF
import NPK.Bar as Bar
from astropy.table import Table 

from scipy.spatial import KDTree 
import scipy.signal as SG
from scipy.interpolate import interp1d



import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
reload(FF)
reload(Extraction)
reload(Wavelength)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''MakeFlat.py

            Creates a flat field image (default name is flat.fits)

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('toflat', type=str, help='Fine correction path')
    parser.add_argument('--toflat2', type=str, help='Fine correction path', default=None)
    parser.add_argument('--max', type=float, help='Max value to divide by', default=None)
    parser.add_argument('--min', type=float, help='Min value to subtract off', default=None)
    parser.add_argument('--maxfrac', type=float, help='Maximum fractional value', default=0.99)
    parser.add_argument('--fillbelow', type=float, help='Values below fillbelow are set to 0', default=0.01)
    parser.add_argument('--outfile', type=str, help='Output filename', default="flat.fits")

    args = parser.parse_args()

    HDU = pf.open(args.toflat)
    if args.toflat2 is not None:
        dat2 = pf.open(args.toflat2)[0].data
        HDU[0].data += dat2*4
        

    min = args.min
    if min is None: min = np.min(HDU[0].data)

    max = args.max
    if max is None: max = np.max(HDU[0].data) * args.maxfrac

    HDU[0].data -= min
    HDU[0].data /= max

    tofill = HDU[0].data < args.fillbelow
    if tofill.any():
        HDU[0].data[tofill] = 0

    HDU.writeto(args.outfile)

