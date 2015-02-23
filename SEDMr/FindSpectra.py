'''

 Finds SEDM spectra based on a sextractor segmentation map

    [c] 2014 Nick Konidaris

'''

import argparse
from collections import namedtuple
import numpy as np
import pylab as pl
import pyfits as pf
from astropy.modeling import models, fitting
from multiprocessing import Pool
import pdb, sys
import NPK.Bar as Bar

from stsci.tools.gfit import gfit1d 



Coord = namedtuple("XY", "x y")

def spanrange(locs):
    ''' Returns a pair of coordinates where index values are true for a 2d array.

    Args:
        locs[d1,d2]: Boolean array of coordinates
    
    Returns:
        [[min x, min y], [max x, max y]]'''

    global XX, YY, xi, yi


    # To avoid expensive creation
    try: XX
    except: 
        xi = np.arange(locs.shape[0])
        yi = np.arange(locs.shape[1])
        XX,YY = np.meshgrid(xi,yi)


    mnx = np.min(XX[locs])
    mny = np.min(YY[locs])
    mxx = np.max(XX[locs])
    mxy = np.max(YY[locs])


    return [Coord(mnx,mny), Coord(mxx,mxy)]



def load_data(segmapfn, objfn):
    '''Returns HDUs from files segmapfn and objfn. load_data is a thin wrapper
        around pyfits'''

    try:
        segmap = pf.open(segmapfn)
    except Exception, e:
        print "Could not open %s: %s" % (segmapfn, e)

    try:
        obj = pf.open(objfn)
    except Exception, e:
        print "Could not open %s: %s" % (segmapfn, e)

    return (segmap, obj)


def find_segments_helper(seg_cnt):
    # Global is for inter process communication
    global segdat, objdat, polyorder
    SegTrace = namedtuple("Trace", "segnum xs ys poly")
    PAD = 2

    the_seg = (segdat == seg_cnt)

    span = spanrange(the_seg)
    mnsr = np.max((span[0].y - PAD, 0))
    mxsr = np.min((span[1].y + PAD, segdat.shape[1]-1))
    y_slc = slice(mnsr, mxsr)
    y_off = (span[0].y+span[1].y)/2.0

    n_el = span[1].x - span[0].x

    if n_el < 50: 
        tr = {"seg_cnt": seg_cnt, "xs": np.array(np.nan), 
            "mean_ys": np.array(np.nan),
            "coeff_ys": np.array(np.nan),
            "trace_sigma": np.nan,
            "ok": False}

        return tr

    means = np.zeros(n_el)
    amps = np.zeros(n_el)
    sds = np.zeros(n_el)
    trace_profile = np.zeros(mxsr-mnsr)

    for i in xrange(n_el):
        XX = i+span[0].x
        profile = np.median(objdat[y_slc, XX-3:XX+3], 1)
        profile -= np.min(profile)

        trace_profile += profile

        xs = np.arange(len(profile)) + span[0].y

        means[i] = np.sum(xs*profile)/np.sum(profile)-PAD
    
    xs = np.arange(n_el) + span[0].x
    poly = np.polyfit(xs, means, polyorder)
    tracefit = gfit1d(trace_profile, par=[0, len(trace_profile)/2., 1.7], quiet=1)

    tr = {"seg_cnt": seg_cnt, "xs": np.array(xs), "mean_ys": np.array(means),
        "coeff_ys": np.array(poly), "ok": True, "trace_sigma": np.abs(tracefit.params[2])}

    print '%4.4i: fwhm=%3.2f pix' % (seg_cnt, np.abs(tracefit.params[2])*2.355)
    return tr

    if plot:
        pl.plot(xs, means, 'x')
        pl.plot(xs, ff(xs))


def find_segments(segmap=None, obj=None, plot=False, order=2):
    '''Find the segments in obj by tracing ridgelines identified in the segmentation map.

    Args:
        segmap[int,int]: Segmentation map image, first segment is 1 .. max(segmap)
        obj[float,float]: Image data to trace over
        plot: Make an example plot (for debugging)
        polyorder: The order of the polynomial used in coeff_ys

    Returns:
        [ (segmap dictionary) ] - Returns a list f with length equal to the max(segmap). segmentation map Dictionaries containing
    {   "seg_cnt": Segment ID number, 
        "xs": List of x positions of trace,
        "mean_ys": Measured Y position (average) of the trace, 
        "coeff_ys": polyfit coefficients to the mean_ys,
        "trace_sigma": Width of trace in pixels (1 sigma),
        "ok": Trace has more than 50 pixels}
        
    '''
    global segdat, objdat, polyorder

    polyorder = order
        
    segdat = segmap[0].data
    objdat = obj[0].data

    PAD = 2

    if plot:
        pl.figure(1)
        pl.clf()
        pl.imshow(objdat, vmin=0,vmax=2000)
        segrange = [500]
    else:
        # First segment begins at 1
        segrange = xrange(1, max(segdat.flatten()))


    p = Pool()
    traces = p.map(find_segments_helper, segrange)
    p.close()

    return traces

def ds9_line_str(startx,starty, endx,endy):
    '''Returns a string representing a line in ds9 '''

    return "line(%s,%s,%s,%s) # line=0 0\n" % (startx,starty,endx,endy)

def write_reports(segmap, obj, segments, outname):
    

    ds9 = "global color=blue\nphysical\n"
    for seg in segments:
        if seg['ok'] == False: continue


        ff = np.poly1d(seg['coeff_ys'])
        for dx in [0,100,200]:
            X1 = seg['xs'][0] + dx
            Y1 = ff(X1)+1
            X2 = X1 + 100
            Y2 = ff(X2)+1

            if np.isfinite(X1) and np.isfinite(X2) and np.isfinite(Y1) and \
                np.isfinite(Y2):
                ds9 += ds9_line_str(X1,Y1, X2,Y2)
    
    fname = "ds9_%s.reg" % outname
    f = open(fname, 'w')
    f.write(ds9)
    f.close()

def write_segments(segmap, obj, segments, outname):
    np.save(objfn + "_segments", segments)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
    '''
    FindSegments.py
    ''')

    parser.add_argument("segmap_fits", type=str, help="Segmentation map .fits file name")
    parser.add_argument("objfn", type=str, help="Dome flat")
    parser.add_argument("outname", type=str, help="Output file name prefix")
    parser.add_argument("--order", type=int, help="Polynomial order", 
                        default=2)
    
    args = parser.parse_args()
    segmapfn = args.segmap_fits
    objfn = args.objfn
    outname = args.outname
    order = args.order

    segmap, obj = load_data(segmapfn, objfn)
    segments = find_segments(segmap, obj, order=order)
    write_segments(segmap, obj, segments, outname)
    write_reports(segmap, obj, segments, outname)


