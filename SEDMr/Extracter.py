
import argparse, os, pdb, sys
import numpy as np
import pylab as pl
import pyfits as pf

from pyraf import iraf

from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d
import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
import NPK.Standards as Stds
reload(Wavelength)
reload(Extraction)


def identify_spectra(spectra, outname=None, low=-np.inf, hi=np.inf, plot=False):


    ms = []
    ixs = []
    segids = []

    ds9 = "physical\n"
    for ix,spectrum in enumerate(spectra):
        if spectrum.__dict__.has_key('spec') and spectrum.spec is not None \
            and spectrum.lamcoeff is not None:
            l,s = spectrum.get_flambda()
            ms.append(np.median(s))
            ixs.append(ix)
            segids.append(spectrum.seg_id)

            X = np.mean(spectrum.xrange)
            Y = np.poly1d(spectrum.poly)(X)

            ds9 += 'point(%s,%s) # point=cross text={%s:%i}\n' % \
                (X,Y, segids[-1],ms[-1])

    ixs = np.array(ixs)
    ms = np.array(ms)
    bgd = np.median(ms)
    sd = np.std(ms)
    segids = np.array(segids)
        
    ms -= bgd

    ok = (ms <= sd*low) | (ms >= sd*hi)
    pl.figure(1)
    pl.clf()
    bgd = (ms > sd*low) & (ms < sd*hi)
    pl.plot(segids[bgd], ms[bgd],'.')
    pl.plot(segids[ok], ms[ok],'x')
    pl.axhline(sd*low,color='orange')
    pl.axhline(sd*hi,color='red')
    pl.xlabel("Sextractor Segment ID number")
    pl.ylabel("Median spectral irradiance [photon/10 m/nm]")
    pl.legend(["Bgd spectra", "Selected spectra"])
    if outname is not None:
        pl.savefig("selected_%s.pdf" % outname)
    if plot:
        pl.show()

    f = open(outname + ".reg", "w")
    f.write(ds9)
    f.close()


    return ixs[ok]
        

def interp_spectra(all_spectra, six, outname=None, plot=False):
    
    l_grid = None
    s_grid = None
    spectra = []
    for ix,spectrum in enumerate(all_spectra):
        if ix not in six: continue

        l,s = spectrum.get_flambda()
        if np.median(s) < 0: pon = -1.0
        else: pon = 1.0

        if l_grid is None:
            l_grid = l
            s_grid = [s*pon]
            lamcoeff = spectrum.lamcoeff
        else:
            fun = interp1d(l,s*pon, bounds_error=False,fill_value=0)
            s_grid.append(fun(l_grid))

        spectra.append(s)
            
            
    medspec = np.mean(s_grid, 0)

    pl.figure(3)
    pl.clf()
    pl.step(l_grid,medspec)
    yl = pl.ylim()
    pl.xlabel('Wavelength [nm]')
    pl.ylabel(r'Spectral irradiance[photon/10 m/nm]')
    if outname is not None: pl.savefig("spec_%s" % outname)
    if plot: pl.show()

    pl.figure(2)
    pl.clf()
    s_grid = np.array(s_grid)
    pl.imshow(s_grid,vmin=yl[0], vmax=yl[1])
    pl.xlabel('Wavelength bin [pixel]')
    pl.colorbar()
    if outname is not None: pl.savefig("allspec_%s" % outname)
    if plot:pl.show()

    
    doc = '''Result contains:
        nm [N float]: Wavelength solution
        ph_10m_nm [N float]: Spectral irradiance of source in units of photon / 10 minute / nm
        spectra [? x K float]: List of all the spectra that participated in
            the formation of ph_10m_nm. By interpolating these objects onto
            a ph_10m_nm and taking the mean, you produce ph_10m_nm
        coefficients [3-5 element float]: Chebyshev coefficents that produce
            nm. Can be evaluated with numpy chebval().
        doc: This doc string
        '''
    result = [{"nm": l_grid, "ph_10m_nm": medspec, "spectra": spectra,
        "coefficients": lamcoeff, 
        "doc": doc}]
    return result

def load_corr():
    corr = pf.open("CORR.npy")

    

def imarith(operand1, op, operand2, result):
    from pyraf import iraf
    iraf.images()

    pars = iraf.imarith.getParList()
    iraf.imcombine.unlearn()

    print "%s %s %s -> %s" % (operand1, op, operand2, result)
    iraf.imarith(operand1=operand1, op=op, operand2=operand2, result=result)

    iraf.imarith.setParList(pars)   
    
def subtract(A,B, outname):
    if os.path.exists(outname):
        return pf.open(outname)

    imarith(A, "-", B, outname)

    return pf.open(outname)

def bgd_level(extractions):
    '''Remove background from extractions'''

    levels = []
    for spectrum in extractions:
        if spectrum.__dict__.has_key('spec') and spectrum.spec is not None \
            and spectrum.lamcoeff is not None:
        
            l, Fl = spectrum.get_flambda()

            levels.append(np.median(Fl))

    bgd = np.median(levels)
    sd = np.std(levels)
    pl.plot(levels,'x')
    pl.axhline(bgd)
    pl.axhline(bgd+sd)
    pl.axhline(bgd-sd)
    pl.ylim(-20*sd-bgd,20*sd-bgd)
    pl.show()

def handle_A(A, fine, outname=None, nsighi=2, standard=None):
    '''Processes A files.

    1. Extracts the spectra from the A file
    2. Identifies the object as some multiple of the sigma'''


    fine = np.load(fine)
    if outname is None:
        outname = "%s" % (A)

    spec = pf.open(A)

    if not os.path.exists(outname + ".npy"):
        E = Wavelength.wavelength_extract(spec, fine, filename=outname)
        np.save(outname, E)
    else:
        E = np.load(outname + ".npy")

    six = identify_spectra(E, hi=nsighi, outname=outname+".pdf")
    res = interp_spectra(E, six, outname=outname+".pdf")
    
    if standard is not None:
        wav = standard[:,0]/10.0
        flux = standard[:,1]

        fun = interp1d(wav, flux, bounds_error=False, fill_value = np.nan)
        correction = fun(res[0]['nm'])/res[0]['ph_10m_nm']

        res[0]['std-correction'] = correction


    np.save("spectrum_" + outname, res)


def handle_corr(filename, outname='CORR.npy', objname=None) :

    if outname is None: outname = 'CORR.npy'
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



def handle_AB(A, B, fine, outname=None, nsiglo=-1, nsighi=1):
    '''Processes A-B files.

    1. Subtracts the two files and creates a A-B
    2. Extracts the spectra from the A-B file
    3. Identifies the object as some multiple of the sigma'''


    fine = np.load(fine)
    if outname is None:
        outname = "%sm%s" % (A,B)

    if not outname.endswith(".fits"): outname = outname + ".fits"
    diff = subtract(A,B, outname)


    if not os.path.exists(outname + ".npy"):
        E = Wavelength.wavelength_extract(diff, fine, filename=outname)
        np.save("extracted_" + outname, E)
    else:
        E = np.load(outname + ".npy")

    six = identify_spectra(E, low=nsiglo, hi=nsighi, outname=outname+".pdf")
    res = interp_spectra(E, six, outname=outname+".pdf")
    np.save("sp_" + outname, res)


parser = argparse.ArgumentParser(description=\
    '''Extracter.py:

        
    ''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--A', type=str, help='FITS A file')
parser.add_argument('--B', type=str, help='FITS B file')
parser.add_argument('fine', type=str, help='Numpy fine wavelength solution')
parser.add_argument('--outname', type=str, help='Prefix output name')
parser.add_argument('--nsiglo', type=float, help='Number sigma to extract below', default=-2)
parser.add_argument('--nsighi', type=float, help='Number sigma to extract above', default=2)
parser.add_argument('--std', type=str, help='Name of standard')

args = parser.parse_args()

if __name__ == '__main__':
    
    if args.outname is not None:
        args.outname = args.outname.rstrip('.npy')

    if args.fine == 'CORR':
        # Take atmospheric correction out and store in a separate file
        handle_corr(args.A, outname=args.outname, objname=args.std)
    
    elif args.A is not None and args.B is not None:
        handle_AB(args.A, args.B, args.fine, outname=args.outname,
            nsiglo=args.nsiglo, nsighi=args.nsighi)

    elif args.A is not None:
        if args.std is None:
            handle_A(args.A, args.fine, outname=args.outname,
                nsighi=args.nsighi)
        else:
            star = Stds.Standards[args.std]
            handle_A(args.A, args.fine, outname=args.outname,
                nsighi=args.nsighi, standard=star)
            
    else:
        print "I do not understand your intent, you must specify --A, at least"
