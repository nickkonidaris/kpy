
import argparse, os, pdb, sys
import numpy as np
import pylab as pl
import pyfits as pf
import scipy.signal as SG

from pyraf import iraf

from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d
import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
import SEDMr.Spectra as SS
import SEDMr.GUI as GUI
import NPK.Standards as Stds
reload(Wavelength)
reload(Extraction)
reload(GUI)
reload(SS)

class Background(object):
    interp_fun = None
    infile = None
    outfile = None
    exptime = None
    lam_nm = None
    spec = None

    def __init__(self, lam_nm=None, infile=None, spec=None,
        outfile=None, exptime=None):

        self.infile = infile 
        self.outfile = outfile 
        self.exptime = exptime
        self.lam_nm = lam_nm
        self.spec = spec

    def __call__(self, x):
        if self.interp_fun is None:
            self.interp_fun = interp1d(self.lam_nm, self.spec, 
                bounds_error=False)
        return self.interp_fun(x)

def background_subtract(KT, objpos=None, radius=100, 
    minl=350, maxl=1000, n_std=9, n_iter=8, smoothing=400):
    
    if objpos is None:
        exclude = KT.KT.query_ball_point(objpos, radius)
    else:
        exclude = []


    lams = []
    specs = []
    for ix in xrange(len(KT.data)):
        e = KT.data[ix]

        if ix in exclude: continue
        if not e.ok: continue
        if e.lamrms > 1: continue
        if e.xrange[1] - e.xrange[0] < 200: continue
        if e.yrange[1] < 0 : continue
        if e.yrange[0] < 0 : continue
        if not np.isfinite(e.yrange[0]): continue
        if not np.isfinite(e.yrange[1]): continue

        try:l,s = e.get_flambda()
        except: continue
        lams.append(l)
        specs.append(s)
    
    exptime = e.exptime
    all_lams = np.array([lam for sublist in lams for lam in sublist])
    all_spec = np.array([spec for sublist in specs for spec in sublist])

    ix = np.argsort(all_lams)
    l,s = all_lams[ix], all_spec[ix]

    ok = (l > minl) & (l < maxl) & np.isfinite(l) & np.isfinite(s)
    knots = np.arange(minl, maxl,.1)
    boxcar = SG.boxcar(smoothing)/smoothing

    nok = len(s[ok])
    for i in xrange(n_iter):
        smoothed = SG.convolve(s[ok], boxcar, mode='same')
        ff = interp1d(l[ok], smoothed, kind='linear', 
            bounds_error=False)

        res = (s - ff(l))*exptime

        std = np.abs(res / np.sqrt(s*exptime))
        ok = (l > minl) & (l < maxl) & (std < n_std) & (np.isfinite(l)) 

        print i, nok, len(s[ok])
        if (float(nok)/len(s[ok]) - 1) < .001: break
        nok = len(s[ok])
        
    n_knots = len(s)/smoothing
    knots = np.arange(minl, maxl, float(maxl-minl)/n_knots)
    bgd = Background(lam_nm=knots, spec=ff(knots), exptime=exptime)

    return bgd

def plot_residuals(KT, bgd, outfile):

    pl.figure(1)
    lg = np.arange(350, 1000)
    residuals = np.zeros((len(lg), len(KT.data)))
    for i in xrange(0,len(KT.data)):
        e = KT.data[i]
        try: l,s = e.get_flambda()
        except: continue
        f = interp1d(l, s, bounds_error=False)
        residuals[:, i] = f(lg) - bgd(lg)

    pl.imshow(residuals, vmin=-.02, vmax=.02)
    pl.colorbar()
    pl.savefig(outfile + ".pdf")

    
    pl.clf()
    pl.plot(l, bgd(l))

    pl.savefig(outfile + "1d.pdf")

def background_subtract_file(infile, outfile, radius=100, minl=350, maxl=1000,
    n_std=9, n_iter=8, smoothing=400):
    exts = np.load(infile)

    KT = SS.Spectra(exts)
    g = GUI.PositionPicker(KT)
    pos = g.picked

    bgd = background_subtract(KT, objpos=pos, radius=radius, minl=minl, 
        maxl=maxl, n_std=9, n_iter=8, smoothing=400)
    bgd.infile = infile
    bgd.outfile = outfile

    np.save(outfile, [bgd])

    plot_residuals(KT, bgd, outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Background.py:

            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('extraction', type=str, help='Numpy extracted file')
    parser.add_argument('outfile', type=str, help='Numpy Background class')
    parser.add_argument('--radius', type=int, default=100, help='Radius of object extraction')
    parser.add_argument('--n_std', type=float, default=9, help='Number of deviations to throw out')
    parser.add_argument('--n_iter', type=int, default=8, help='Number of cleaning iterations')
    parser.add_argument('--smoothing', type=int, default=400, help='Amount of smoothing of background')

    args = parser.parse_args()

    
    background_subtract_file(args.extraction, args.outfile, radius=args.radius,
        n_std=args.n_std, n_iter=args.n_iter, smoothing=args.smoothing)


