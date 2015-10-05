
import argparse
import numpy as np
import pylab as pl
import pyfits as pf
from scipy.interpolate import interp1d
import sys

outdir = "/scr2/npk/sedm/reduced/"
VER = "2015A"

 
def archive(objname):

    corrname = "std-correction.npy"
    if not os.path.isfile(corrname):
        corrname= '/scr2/npk/sedm/OUTPUT/2015mar25/std-correction.npy'

    ss = np.load(specname)[0]
    corr = np.load(corrname)[0]

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
    pl.ion()
    pl.show()


    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Check.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('name', type=str, help='Object name to archive')


    archive(name)


