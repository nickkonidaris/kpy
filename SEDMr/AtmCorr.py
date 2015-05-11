
import argparse, os, sys
import numpy as np
import pylab as pl
import pyfits as pf

from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d

def handle_summary(outname=None, filelist=[]):
    if outname is None: outname = 'correction.npy'

    keepers = []
    for file in filelist: 
        print file
        f = np.load(file)[0]
        if np.nanmin(f['std-correction']) < 9:
            keepers.append(f)
            print "Keeping %s" % file


    corl = keepers[0]['nm'].copy()
    cor = np.zeros((len(corl), len(keepers)))
    cor[:,0] = keepers[0]['std-correction'].copy()
    for ix, keeper in enumerate(keepers[1:]):
        f = interp1d(keeper['nm'], keeper['std-correction'], bounds_error=False, 
            fill_value = np.nan)

        cor[:,ix] = f(corl)

    cs = np.nanmean(cor,1)
    ccs = chebfit(corl, cs, 6)

    cor = [{"nm": corl, "cor": cs, "coeff": ccs}]
    np.save(outname, cor)

def handle_corr(filename, outname='std-correction.npy', objname=None) :

    if outname is None: outname = "corr_" + filename
    dat = np.load("spectrum_" + filename)[0]
    
    pl.figure(1)
    pl.clf()
    pl.xlim(365, 1000)
    pl.semilogy(dat['nm'], dat['std-correction'])
    pl.ylim(1,1e4)
    pl.grid(True)
    pl.savefig("corr_" + filename.rstrip(".npy") + ".pdf")

    if os.path.exists(outname):
        res = np.load(outname)
    else:
        res = []
    
    res.append({objname: {'nm': dat['nm'], 'cor': dat['std-correction']}})

    np.save(outname, res)




parser = argparse.ArgumentParser(description=\
    '''

        
    ''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('process', type=str, help='Process [CORR|SUM]')
parser.add_argument('--A', type=str, help='FITS A file')
parser.add_argument('--outname', type=str, help='Prefix output name')
parser.add_argument('--std', type=str, help='Name of standard')
parser.add_argument('--files', type=str, metavar='file', nargs='+',
    help='Name of standard')


args = parser.parse_args()

if __name__ == '__main__':
    
    if args.process == 'CORR':
        # Take atmospheric correction out and store in a separate file
        handle_corr(args.A, outname=args.outname, objname=args.std)

    if args.process == 'SUM':
        handle_summary(outname=args.outname, filelist=args.files)
 




