
import argparse, os, pdb, sys
import numpy as np
import pylab as pl
import pyfits as pf
import scipy.signal as SG
from scipy.spatial import KDTree 

from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d
import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
import SEDMr.Spectra as SS

import sys



class FitsCube(object):
    
    # Following are numpy arrays
    cube = None
    rect_cube = None

    cube_snr = None
    rect_cube_snr = None

    locations = None
    header = None

    
    def load(fname):
        FF = pf.open(fname)

        header = FF[0].header
        cube = FF[0].data
        cube_snr = FF[1].data
        locations = FF[2].data
        rect_cube = FF[3].data
        rect_cube_snr = FF[4].data


        
        

        


