
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


def measure_flexure_x(fine, HDUlist, plot=True, dY=0,
    skyline=589.0, lamstart=1000.0, lamratio=239./240., lamlen=250,
    extract_width=3, skywidth=9, outfile='dX'):
    '''Measures flexure in X direction, returns pixel offset

    Args:
        fine: List of Extraction object, the fine loc + wave solution for
            each spectrum
        HDUlist: Pyfits object for the spectrum to measure
        plot: Plot + save results to a file
        dY: the measured pixel flexure in Y direction to account for
        skyline(float): The night skyline to centroid on in nm
        skywidth(float): Fit gaussian to the ROI of (skyline-skywidth to
            skyline+skywidth) in nm.
        extract_width(int): Number of pixels to extract spectrum around
        - See Wavelength.fiducial spectrum for following:
        lamstart: Wavelength to start the grid on, default 1000 nm
        lamratio: Resolution of sed machine
        lamlen: Length of spectrum

    Returns:
        Offset number of pixels in X direction.
    '''
    
    dat = HDUlist[0].data
    exptime = HDUlist[0].header['EXPTIME']

    spec_ixs = np.arange(500, 1200, 10)
    lamgrid = Wavelength.fiducial_spectrum(lamstart=lamstart,
        lamratio=lamratio, len=lamlen)

    specgrid = np.zeros((len(lamgrid), len(spec_ixs)))
    
    for i,ix in enumerate(spec_ixs):
        f = fine[ix]

        if not f.ok: continue
        if f.lamrms > 1: continue
        if f.xrange[1] - f.xrange[0] < 200: continue

        spec = np.zeros(f.xrange[1] - f.xrange[0])
        yfun = np.poly1d(f.poly)

        for jx,xpos in enumerate(np.arange(f.xrange[0], f.xrange[1])):
            ypos = yfun(xpos)

            try:spec[jx] = np.sum(dat[ypos-extract_width:ypos+extract_width,
                    xpos])
            except: continue

        try:ll = f.get_lambda_nm()
        except: continue
        specfun = interp1d(ll, spec, bounds_error=False)
        specgrid[:,i] = specfun(lamgrid)
            
    skyspec = np.median(specgrid, axis=1)
    pl.step(lamgrid, skyspec, where='mid')

    roi = (lamgrid>skyline-skywidth) & (lamgrid<skyline+skywidth)
    ffun = FF.mpfit_residuals(FF.gaussian4)
    parinfo= [
        {'value': np.max(skyspec[roi]), 'limited': [1,0], 
            'limits': [0, 0]},
        {'value': skyline}, 
        {'value': 3}, 
        {'value': np.min(skyspec[roi]), 'limited': [1,0],
            'limits': [0,0]}]
    fit = FF.mpfit_do(ffun, lamgrid[roi], skyspec[roi], parinfo)
    pl.plot(lamgrid, FF.gaussian4(fit.params, lamgrid))

    pl.savefig(outfile + ".pdf")

    dXnm = fit.params[1] - skyline

    print "dX = %3.2f nm shift" % dXnm


    return dXnm

 
def measure_flexure_y(fine, HDUlist, profwidth=5, plot=False, outname=None):
    
    dat = HDUlist[0].data
    exptime = HDUlist[0].header['EXPTIME']

    profs = []
    xx = np.arange(profwidth*2)

    for ix in np.arange(500, 1200, 10):
        f = fine[ix]
        profile = np.zeros(profwidth*2)

        if not f.ok: continue
        if f.lamrms > 1: continue
        if f.xrange[1] - f.xrange[0] < 200: continue

        yfun = np.poly1d(f.poly)
        for xpos in np.arange(f.xrange[0], f.xrange[1]):
            try:    ypos = int(np.round(yfun(xpos)))
            except: continue
            try: profile += dat[ypos-profwidth:ypos+profwidth, xpos]
            except: continue
            
        profs.append(profile - np.min(profile))

    if plot: pl.figure(1)
    ffun = FF.mpfit_residuals(FF.gaussian4)
    parinfo= [{'value': 1}, {'value': profwidth}, {'value': 2}, 
        {'value': 0}, {'value': 0}]

    profposys = []
    for prof in profs:
        if plot: pl.step(xx, prof)
        parinfo[0]['value'] = np.max(prof)
        fit = FF.mpfit_do(ffun, xx, prof, parinfo)
        if plot: pl.plot(xx, FF.gaussian4(fit.params, xx))
        profposys.append(fit.params[1] - profwidth-1)
    if plot: pl.show()
    profposys = np.array(profposys)

    mn = np.mean(profposys)
    sd = np.std(profposys)
    ok = np.abs(profposys - mn)/sd < 3
    required_shift = np.mean(profposys[ok])
    print "dY = %3.2f pixel shift" % required_shift
        
    return required_shift



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Flexure.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('fine', type=str, help='Fine correction path')
    parser.add_argument('infile', type=str, help='Path to FITS file to refit')
    parser.add_argument('--profwidth', type=float, help='Profile width to extract for Y flexure', default=5)
    parser.add_argument('--skyline', type=float, help='skyline position in nm to measure X flexure', default=589.0)
    parser.add_argument('--lamstart', type=float, help='Wavelength to start interpolating grid ', default=1000.0)
    parser.add_argument('--lamratio', type=float, help='Wavelength resolution for interpolating grid', default=239.0/240.0)
    parser.add_argument('--lamlen', type=int, help='Wavelength grid length', default=250)
    parser.add_argument('--extract_width', type=int, help='Extraction width for spectrum in the Y direction (for X flexure)', default=3)
    parser.add_argument('--skywidth', type=float, help='Wavelength range to search over for X flexure measurement. Range is (skyline-sky_width : skyline+sky_width', default=9)
    parser.add_argument('--outfile', type=str, help='Output filename', default="flexure.npy")

    args = parser.parse_args()

    fine = np.load(args.fine)
    HDU = pf.open(args.infile)
    dy = measure_flexure_y(fine, HDU, profwidth=args.profwidth)
    dx = measure_flexure_x(fine, HDU, dY=dy,
        skyline=args.skyline,
        lamstart=args.lamstart,
        lamratio=args.lamratio,
        lamlen=args.lamlen,
        extract_width=args.extract_width,
        skywidth=args.skywidth,
        outfile=args.outfile)

    res = [{'fine_name': args.fine,
        'infile_name': args.infile,
        'profwidth': args.profwidth,
        'skyline': args.skyline,
        'lamstart': args.lamstart,
        'lamratio': args.lamratio,
        'lamlen': args.lamlen,
        'extract_width': args.extract_width,
        'skywidth': args.skywidth,
        'dXnm': dx,
        'dYpix': dy}]

    np.save(args.outfile, res)
