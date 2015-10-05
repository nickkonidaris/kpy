
'''
    25 June 2015- I'm testing multi-peak Lorentz fitting for extracting spectra
'''


import argparse, copy, os, pdb, sys
import numpy as np
import pylab as pl
import pyfits as pf
import scipy.signal as SG
import itertools
import collections

from scipy.optimize import leastsq
from scipy.spatial import KDTree 


from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d
import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
import SEDMr.Spectra as SS
import SEDMr.GUI as GUI
import NPK.Standards as Stds


def lorentz(p,x):
    '''
    offset y0  = p[0]
    area A   = p[1]
    w   = p[2]
    xc  = p[3]'''
    y0  = p[0]
    A   = np.abs(p[1])
    w   = np.abs(p[2])
    xc  = p[3]
        
    return y0 + 2*A/np.pi * w / ((x - xc)**2 + w**2)


def lorentz_n(p, xcs, x):

    y0  = p[0]
    w   = np.abs(p[1])
    As  = p[2:]

    if len(As) != len(xcs):
        import pdb
        pdb.set_trace()
        raise Exception("Mismatch in Lorentz parameter length")

    ret = y0
    for i in xrange(len(As)):
        ret += 2*As[i]/np.pi * w/((x-xcs[i])**2 + w**2)

    return ret

def errfunc(p, xcs, x, z):
    '''


    '''

    #### INLINE lorentz_n function BY HAND
    y0  = p[0]
    w   = np.abs(p[1])
    As  = p[2:]

    if len(As) != len(xcs):
        import pdb
        pdb.set_trace()
        raise Exception("Mismatch in Lorentz parameter length")

    ret = y0
    for i in xrange(len(As)):
        ret += 2*As[i]/np.pi * w/((x-xcs[i])**2 + w**2)
    #### END INLINE

    return  ret-z
        
def handle():
    
    fine = np.load('fine.npy')
    traces = np.load('dome.fits_segments.npy')
    FF = pf.open('STD-G193-74-airmass1.1_v1.fits')
    FF = pf.open('b_ifu20150325_13_18_58.fits')

    header = FF[0].header
    data = FF[0].data

    FF2 = pf.open("seg_dome.fits")
    segs = FF2[0].data

    seg_id = 1817

    sids = []
    xs = []
    ys = []
    good_segs = []
    Results = fine[:]

    for seg in traces:
        if not seg['ok']: continue
        y = np.nanmean(seg['mean_ys'])
        if y != y: continue
        
        sids.append(seg['seg_cnt'])
        good_segs.append(seg)
        xs.append(np.nanmean(seg['xs']))
        ys.append(y)



    for res in Results:
        res.nmeas = np.zeros(265+80)
        res.meas = np.zeros((265+80,10))
        res.widths = np.zeros((265+80, 10))
        res.y0s = np.zeros((265+80, 10))

    sids = np.array(sids)
    dat = np.array([xs,ys],dtype=np.float).T
    KT = KDTree(dat)

    seg_ix = np.where(sids==seg_id)[0][0]
    seg = good_segs[seg_ix]

    y = np.nanmean(seg['mean_ys'])
    print y

    HALFWIDTH = 30

    yvals = np.arange(np.round(y-HALFWIDTH-1), np.round(y+HALFWIDTH+1))

    #for x in np.arange(*seg.xrange):
    results = []
    ResImage = np.zeros((2048, 2048))

    for x_cnt, x in enumerate(seg['xs']):
        spec = data[yvals.astype(np.int), x]
        to_search = KT.query_ball_point((x,y), 250)

        seg_ids = []
        locations = []
        print

        for candidate_ix in to_search:
            candidate = good_segs[candidate_ix]
            xr = candidate['xs']
            x_location = np.argmin(np.abs(xr - x))

            ff = np.poly1d(candidate['coeff_ys'])

            if (xr[0]-80) <= x <= (xr[0] + 265) and np.abs(y - ff(x_location)) < (HALFWIDTH-1):
                seg_ids.append(candidate['seg_cnt'])
                locations.append(candidate['mean_ys'][x_location])
                print "**: ", candidate['seg_cnt'], locations[-1]
        
        intf = interp1d(yvals, spec, bounds_error=0, fill_value=0)
        y0 = np.percentile(spec, 5)
        w = 1
        As  = intf(locations)
        xcs = np.copy(locations)

        ps = np.hstack([y0, w, As])
        
        fitps = leastsq(errfunc, ps, args=(xcs, yvals, spec))
        fitps = fitps[0]

        ResImage[np.round(yvals).astype(np.int), x] = lorentz_n(fitps, xcs, np.round(yvals))

        fit = []
        for ix in xrange(len(seg_ids)):
            r = Results[seg_ids[ix]-1]

            into_x = x - r.xrange[0] 
            r.meas[into_x, r.nmeas[into_x]] = fitps[2+ix]
            r.widths[into_x, r.nmeas[into_x]] = fitps[1]
            r.y0s[into_x, r.nmeas[into_x]] = fitps[0]
            r.nmeas[into_x] += 1

    OUT = pf.PrimaryHDU(ResImage)
    OUT.writeto("DIMG.fits")
    import pdb
    pdb.set_trace()

    ext_1817 = Results[1816].meas[:,0]

    np.save("ext_1817.npy", ext_1817)

if __name__ == '__main__':
    print '''
    cd /scr2/npk/sedm/new_process/2015mar25
    pharos% ~/spy /scr2/npk/PYTHON/SEDMr/ExtracterTest25jun2015.py
    '''


    handle()
