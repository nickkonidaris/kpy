
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


from numpy.polynomial.chebyshev import chebfit, chebval

import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
reload(FF)
reload(Extraction)
reload(Wavelength)


def measure_flat(fine, HDUlist, 
        lamstart=1000,
        lamratio=239./240.,
        lamlen=250,
        extract_width=3,
        order=7,
        outfile='flat.npy'):

    dat = HDUlist[0].data

    lamgrid = Wavelength.fiducial_spectrum(lamstart=lamstart,
        lamratio=lamratio, len=lamlen)

    specgrid = np.zeros((len(lamgrid), len(fine)))
    
    for i,f in enumerate(fine):

        if not f.ok: continue
        if f.lamrms > 1: continue
        if f.xrange[1] - f.xrange[0] < 200: continue

        spec = np.zeros(f.xrange[1] - f.xrange[0])
        yfun = np.poly1d(f.poly)

        for jx,xpos in enumerate(np.arange(f.xrange[0], f.xrange[1])):
            ypos = np.round(yfun(xpos))

            try:spec[jx] = np.sum(dat[ypos-extract_width:ypos+extract_width,
                    xpos])
            except: continue

        try:ll = chebval(np.arange(len(spec)), f.lamcoeff)
        except: continue
        specfun = interp1d(ll, spec, bounds_error=False)
        specgrid[:,i] = specfun(lamgrid)
            
    flatspec = np.median(specgrid, axis=1)

    chebs = np.zeros((order+1, len(fine)))
    for i in np.arange(len(fine)):
        chebs[:,i] = chebfit(lamgrid, specgrid[:,i] / flatspec, order)

    return chebs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Flexure.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('fine', type=str, help='Fine correction path')
    parser.add_argument('infile', type=str, help='Path to dome flat')
    parser.add_argument('--lamstart', type=float, help='Wavelength resolution for interpolating grid', default=1000.0)
    parser.add_argument('--lamratio', type=float, help='Wavelength resolution for interpolating grid', default=239.0/240.0)
    parser.add_argument('--lamlen', type=int, help='Wavelength grid length', default=250)
    parser.add_argument('--extract_width', type=int, help='Extraction width for spectrum in the Y direction (for X flexure)', default=3)
    parser.add_argument('--order', type=int, help='Chebyshev polynomial order', default=7)
    parser.add_argument('--outfile', type=str, help='Output filename', default="flexure.npy")

    args = parser.parse_args()

    fine = np.load(args.fine)
    HDU = pf.open(args.infile)
    flat = measure_flat(fine, HDU, 
        lamstart=args.lamstart,
        lamratio=args.lamratio,
        lamlen=args.lamlen,
        extract_width=args.extract_width,
        order=args.order,
        outfile=args.outfile)

    np.save(args.outfile, flat)
