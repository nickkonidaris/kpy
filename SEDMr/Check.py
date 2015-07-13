
import argparse
import numpy as np
import pylab as pl
import pyfits as pf
from scipy.interpolate import interp1d
import sys

OUTDIR = '/scr2/npk/sedm/reduced/2015A'

 
def checkSpec(specname):

    ss = np.load(specname)[0]

    print "Plotting spectrum in %s" % specname
    print "Extraction radius: %1.2f" % ss['radius_as']
    print ss.keys()

    corr = np.load('atm-corr.npy')[0]
    corf = interp1d(corr['nm'],corr['correction'], bounds_error=False,
        fill_value=1.0)

    ec = 0
    ext = None
    if ss.has_key('extinction_corr'):
        ext = ss['extinction_corr']
        ec = np.median(ext)
    elif ss.has_key('extinction_corr_A'):
        ext = ss['extinction_corr_A']
        ec = np.median(ext)

    et = ss['exptime']
    pl.title("%s\n(airmass corr factor ~ %1.2f Exptime: %i)" % (specname, ec, et))
    pl.xlabel("Wavelength [nm]")
    pl.ylabel("erg/s/cm2/ang")
    
    pl.step(ss['nm'], ss['ph_10m_nm']*corf(ss['nm']), linewidth=2)
    try: pl.step(ss['skynm'], ss['skyph']*corf(ss['skynm']))
    except: pl.step(ss['nm'], ss['skyph']*(ss['N_spaxA']+ss['N_spaxB'])*
        corf(ss['nm']))

    try: pl.step(ss['nm'], np.sqrt(np.abs(ss['var']))*corf(ss['nm']))
    except: pass

    pl.legend(['obj', 'sky', 'err'])
    pl.xlim(360,1100)
    pl.grid(True)
    pl.ioff()
    pl.show()

def checkCube(cubename):
    ''' Plot a datacube for checking'''
    
    cc = np.load(cubename)

    Xs = [c.X_as for c in cc]
    Ys = [c.Y_as for c in cc]
    Ss = [c.trace_sigma for c in cc]

    pl.figure(1)
    pl.scatter(Xs, Ys, marker='H', linewidth=0, s=50, c=Ss, vmin=0.8, vmax= 2)
    pl.title("This should show a regular hexagonal grid of cube positions")
    pl.xlim(-25,25)
    pl.ylim(-25,25)

    pl.colorbar(label='RMS trace width in pixel')
    pl.xlabel("X position [as]")
    pl.ylabel("Y position [as]")
    pl.grid(True)
    pl.ioff()
    pl.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Check.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('name', type=str, help='Name of object to archive')



