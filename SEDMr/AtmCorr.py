
import argparse, os, sys
import numpy as np
import pylab as pl
import pyfits as pf
import datetime
import os
import sets

from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d

import Wavelength
import NPK.Standards as SS

def handle_create(outname=None, filelist=[]):
    ''' Create standard star correction. Units are erg/s/cm2/Ang '''

    if outname is None: outname='atm-corr.npy'

    ll = Wavelength.fiducial_spectrum()
    spectra = []
    corrs = []
    corr_vals =[]
    for file in filelist: 
        ''' Filename is STD-nnnnn_obs* '''

        try: data = np.load(file)[0]
        except:
            raise Exception("Not passed a spectrum file, the file should start with spectrum_")
        
        # Convert file named spectrum_STD-nnnn_obs* to a name compaitable
        # with Standards.py
        pred = file[13:].rstrip()
        pred = pred.split("_")[0]
        pred = pred.lower().replace("+","").replace("-","_")
        print pred, pred in SS.Standards

        if pred not in SS.Standards:
            raise exception("File named '%s' is reduced to '%s' and no such standard seems to exist."  % (file, pred))


        ff = interp1d(data['nm'], data['ph_10m_nm'], bounds_error=False)
        spectra.append(ff(ll))


        # Now process the standard
        standard = SS.Standards[pred]
        wav = standard[:,0]/10.0
        flux = standard[:,1]
        std_ff = interp1d(wav, flux, bounds_error=False, fill_value=np.nan)
        correction = std_ff(ll)/ff(ll)

        ROI = (ll > 600) & (ll < 850)
        corr_vals.append(np.median(correction[ROI]))
        correction /= corr_vals[-1]

        corrs.append(correction)

    corrs = np.array(corrs) 
    erg_s_cm2_ang = corrs * np.median(corr_vals) * 1e-16


    
    roi = (ll > 900) & (ll < 1000)
    the_corr = np.median(erg_s_cm2_ang,0)
    if not np.isfinite(the_corr).all():
        # Fit polynomial to extrapolate correction
        redend = (ll > 800) & (ll < 915)
        ss = np.poly1d(np.polyfit(ll[redend], the_corr[redend], 2))

        redend = (ll>910)
        the_corr[redend] = ss(ll[redend])


    # Now clean over the Balmer series
    balmers = [656.3, 486.1, 434.0, 410.2, 397.0]
    for balmer in balmers:
        eps = balmer * 0.02
        line_ROI = sets.Set(np.where(np.abs((ll-balmer)/balmer) < 0.03)[0])
        broad_ROI = sets.Set(np.where(np.abs((ll-balmer)/balmer) < 0.06)[0])
        around_line_ROI = list(broad_ROI - line_ROI)
        fit = np.poly1d(np.polyfit(ll[around_line_ROI], 
            the_corr[around_line_ROI], 1))
        to_fix = list(line_ROI)
        the_corr[to_fix] = fit(ll[to_fix])



    # Plot data
    pl.figure(1)
    pl.clf()
    pl.grid(True)
    pl.ylim(1e-20,1e-16)
    pl.semilogy(ll, the_corr, linewidth=4)
    for ix,e in enumerate(erg_s_cm2_ang):
        pl.semilogy(ll, e*corr_vals[ix]/np.mean(corr_vals))

    pl.xlabel("Wavelength [nm]")
    pl.ylabel("Correction [erg/s/cm cm/Ang]")
    pl.title("Correct ph/10 m/nm to erg/s/cm2/Ang")
    pl.savefig("Standard_Correction.pdf")
    print np.mean(corr_vals) * 1e-16, np.std(corr_vals)*1e-16

    # Construct result
    res = {"nm": ll,
        "correction": the_corr,
        "doc": "Correct ph/10 m/nm to erg/2/cm2/ang",
        "Nspec": len(corrs),
        "correction_std": np.nanstd(erg_s_cm2_ang,0),
        "outname": outname,
        "files": filelist,
        "when": '%s' % datetime.datetime.now(),
        "user": os.getlogin()
        }

    np.save(outname, [res])
    return res

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




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('process', type=str, help='Process [CORR|SUM|CREATE]')
    parser.add_argument('--A', type=str, help='FITS A file')
    parser.add_argument('--outname', type=str, help='Prefix output name')
    parser.add_argument('--std', type=str, help='Name of standard')
    parser.add_argument('--files', type=str, metavar='file', nargs='+',
        help='Name of standard')


    args = parser.parse_args()


    
    if args.process == 'CORR':
        # Take atmospheric correction out and store in a separate file
        handle_corr(args.A, outname=args.outname, objname=args.std)

    if args.process == 'SUM':
        handle_summary(outname=args.outname, filelist=args.files)

    if args.process == 'CREATE':
        # Create the atmospheric correction.
        handle_create(outname=args.outname, filelist=args.files)
 




