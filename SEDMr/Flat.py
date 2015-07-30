
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
from SEDMr.FlatCorrection import FlatCorrection as FC

reload(FF)
reload(Extraction)
reload(Wavelength)


def measure_flat(extraction, meta, 
        lamstart=700,
        lamend=850,
        outfile='flat.npy'):


    corrections = []
    Xs = []
    Ys = []
    for i,e in enumerate(extraction):
        fc = FC(seg_id=e.seg_id)
        corrections.append(fc)

        if not e.ok: continue

        try: l,f = e.get_flambda()
        except: continue

        X = np.argmin(np.abs(l-656))
        Xs.append(e.xrange[0] + X)
        Ys.append(np.mean(e.yrange))

        ROI = (l>lamstart) & (l <= lamend)
        fc.correction = np.nanmean(f[ROI])

    vals = [f.get_correction(0) for f in corrections]
    medval = np.median(vals)

    Ss = []
    for c in corrections: 
        try: 
            c.correction /= medval
            if c.correction > 2:
                c.correction = 1.0
            Ss.append(c.correction)
        except: pass
    
    pl.figure(1)
    pl.clf()
    pl.scatter(Xs, Ys, c=Ss,s=70,linewidth=0, vmin=0.8,vmax=1.2,marker='h')
    pl.xlim(-100,2048+100)
    pl.ylim(-100,2048+100)
    pl.colorbar()
    pl.xlabel("X pixel")
    pl.ylabel("Y pixel")
    pl.title("Single correction from %s to %s from %s" % (lamstart, lamend,
        meta['outname']))
    pl.savefig("flat-field-values.pdf")

    return corrections

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('infile', type=str, help='Path to dome flat')
    parser.add_argument('--lamstart', type=float, help='Wavelength resolution for interpolating grid', default=700.0)
    parser.add_argument('--lamend', type=float, help='Wavelength resolution for interpolating grid', default=900.0)
    parser.add_argument('--outfile', type=str, help='Output filename', default="flat-dome-700to900.npy")

    args = parser.parse_args()

    ext, meta = np.load(args.infile)
    flat = measure_flat(ext, meta,
        lamstart=args.lamstart,
        lamend=args.lamend)

    np.save(args.outfile, flat)
