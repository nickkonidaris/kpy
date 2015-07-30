
import argparse, copy, os, pdb, sys
import numpy as np
import pylab as pl
import pyfits as pf
import scipy.signal as SG
import itertools



from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d
import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
import SEDMr.Spectra as SS
import SEDMr.GUI as GUI
import NPK.Standards as Stds
import NPK.Atmosphere as Atm

def identify_spectra_gui(spectra, outname=None, radius=2):
    ''' Returns index of spectra picked in GUI.

    NOTE: Index is counted against the array, not seg_id'''
    
    print "Looking in a %s as radius" % radius
    pl.ioff()
    KT = SS.Spectra(spectra)
    g = GUI.PositionPicker(KT, bgd_sub=True, radius_as=radius)
    pos  = g.picked

    kix = KT.KT.query_ball_point(pos, radius)

    pl.close()

    return KT.good_positions[kix], pos

def identify_bgd_spectra(spectra, pos, inner=3, outer=6):
    KT = SS.Spectra(spectra)

    objs = KT.good_positions[KT.KT.query_ball_point(pos, r=inner)]
    skys = KT.good_positions[KT.KT.query_ball_point(pos, r=outer)].tolist()

    for o in objs:
        if o in skys: skys.remove(o)

    return skys




def identify_spectra(spectra, outname=None, low=-np.inf, hi=np.inf, plot=False):

    raise Exception("This code is outdated")

    ms = []
    ixs = []
    segids = []

    ds9 = "physical\n"
    for ix,spectrum in enumerate(spectra):
        if spectrum.__dict__.has_key('spec') and spectrum.spec is not None \
            and spectrum.lamcoeff is not None:
            ixs.append(ix)
            segids.append(spectrum.seg_id)

            try: l,s = spectrum.get_counts(the_spec='specw')
            except: 
                ms.append(0)
                continue


            ms.append(np.median(s))
            X = spectrum.X_as
            Y = spectrum.Y_as

            ds9 += 'point(%s,%s) # point=cross text={%s:%4.2f}\n' % \
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
    pl.grid(True)
    if outname is not None:
        pl.savefig("selected_%s.pdf" % outname)
    if plot:
        pl.show()

    f = open(outname + ".reg", "w")
    f.write(ds9)
    f.close()

    to_image(outname)


    return ixs[ok]

def to_image(spectra, meta, outname, posA=None, posB=None, radius=None):
    ''' Convert spectra list into image_[outname].pdf '''
    def Strength(x):
        if x.xrange is None: return None
        if x.lamcoeff is None: return None
        if x.specw is None: return None
        ix = np.arange(*x.xrange)
        ll = chebval(ix, x.lamcoeff)
        OK = (ll > 450) & (ll < 700)

        if OK.any():
            return np.sum(x.specw[OK])
        else:
            return None

    sig = [Strength(seg) for seg in spectra]

    XS = np.array([seg.X_as for seg in spectra])
    YS = np.array([seg.Y_as for seg in spectra])
    sig = np.array(sig, dtype=np.float)


    pl.clf()
    pl.ylim(-20, 20)
    pl.xlim(-20, 20)
    pl.grid(True)
    if posA is not None:
        pl.axvline(posA[0], color='black', linewidth=.5)
        pl.axhline(posA[1], color='black', linewidth=.5)
    if posB is not None:
        pl.axvline(posB[0], color='black', linewidth=.5)
        pl.axhline(posB[1], color='black', linewidth=.5)
    pl.scatter(XS, YS, c=sig,s=50,marker='H',linewidth=0)

    pl.xlabel("X [as]")
    pl.ylabel("Y [as]")
    pl.title(meta['outname'])
    pl.colorbar()
    pl.savefig("image_%s.pdf" % outname)
    pl.close()

        
def c_to_nm(coefficients, pix, offset=0):
    
    t = coefficients[:]
    t[0] += offset
    return chebval(pix, t)

def interp_spectra(all_spectra, six, sign=1., outname=None, plot=False,
    corrfile=None, dnm=0, onto=None):
    '''Interp spectra onto common grid

    Args:
        all_spectra:
        six:
        dnm: Offset (usually for flexure) in nm'''
    
    l_grid = onto
    s_grid = []
    lamcoeff = None
    #for ix,spectrum in enumerate(all_spectra):
    for ix in six:
        spectrum = all_spectra[ix]

        l,s = spectrum.get_counts(the_spec='specw')
        pix = np.arange(*spectrum.xrange)

        if spectrum.mdn_coeff is not None: cs = spectrum.mdn_coeff
        else: cs = spectrum.lamcoeff
        l = c_to_nm(cs, pix, offset=dnm)
        if l.max() - l.min() < 300: continue

        pon = sign

        if l_grid is None:
            l_grid = l
            s_grid.append(s*pon)
            lamcoeff = spectrum.lamcoeff
        else:
            fun = interp1d(l,s*pon, bounds_error=False,fill_value=0)
            s_grid.append(fun(l_grid))

            
            
    medspec = np.mean(s_grid, 0)


    pl.figure(3)
    pl.clf()
    pl.step(l_grid,medspec)
    yl = pl.ylim()
    pl.xlabel('Wavelength [nm]')
    pl.ylabel(r'Spectral irradiance[photon/10 m/nm]')
    pl.grid(True)
    if outname is not None: pl.savefig("spec_%s" % outname)
    if plot: pl.show()

    pl.figure(2)
    pl.clf()
    s_grid = np.array(s_grid)
    pl.imshow(s_grid,vmin=yl[0], vmax=yl[1])
    pl.xlabel('Wavelength bin [pixel]')
    pl.colorbar()
    pl.grid(True)
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
        corrected-spec [N float]: ph_10m_nm * Atmospheric correction, if 
            available
        doc: This doc string
        '''
    result = [{"nm": l_grid, "ph_10m_nm": medspec, "spectra": s_grid,
        "coefficients": lamcoeff, 
        "doc": doc}]

    CC = None
    if corrfile is not None:
        try: CC = np.load(corrfile)[0]
        except: CC = None

    if CC is not None:
        corrfun = chebval(l_grid, CC['coeff'])
        corrfun /= np.nanmin(corrfun)
        corrfun = interp1d(CC['nm'], CC['cor'], bounds_error=False, fill_value=np.nan)
        corrfun = corrfun(l_grid)
        result[0]['corrected-spec'] = medspec * corrfun
        pl.figure(4)
        pl.clf()
        pl.step(l_grid,medspec*corrfun)
        pl.ylim(yl[0], yl[1]*20)
        pl.xlabel('Wavelength [nm]')
        pl.ylabel(r'Spectral irradiance[photon/10 m/nm] x Atm correction')
        pl.grid(True)
        if outname is not None: pl.savefig("corr_spec_%s" % outname)
        if plot: pl.show()


    pl.figure(2)


    return result

def load_corr():
    corr = pf.open("CORR.npy")

    

def imarith(operand1, op, operand2, result, doAirmass=False):
    from pyraf import iraf
    iraf.images()

    pars = iraf.imarith.getParList()
    iraf.imcombine.unlearn()

    try: os.remove(result)
    except: pass

    print "%s %s %s -> %s" % (operand1, op, operand2, result)
    iraf.imarith(operand1=operand1, op=op, operand2=operand2, result=result)
    iraf.imarith.setParList(pars)   
    if doAirmass:
        # Adjust FITS header
        with pf.open(operand1) as f:
            am1 = f[0].header['airmass']
        with pf.open(operand2) as f:
            am2 = f[0].header['airmass']

        of = pf.open(result)
        of[0].header['airmass1'] = am1
        of[0].header['airmass2'] = am2
        of.writeto(result, clobber=True)

def gunzip(A, B):
    if A.endswith(".gz"):
        os.system("gunzip %s" % A)
    if B.endswith(".gz"):
        os.system("gunzip %s" % B)

    return A.rstrip(".gz"), B.rstrip(".gz")

def gzip(A,B):
    if not A.endswith(".gz"):
        os.system("gzip %s" % A)
    if not B.endswith(".gz"):
        os.system("gzip %s" % B)

    return A+".gz", B+".gz"
    
def add(A,B, outname):
    A,B = gunzip(A,B)
    imarith(A, "+", B, outname)
    gzip(A,B)

    return pf.open(outname)

def subtract(A,B, outname):
    if os.path.exists(outname):
        return pf.open(outname)

    A,B = gunzip(A,B)
    imarith(A, "-", B, outname, doAirmass=True)
    A,B = gzip(A,B)

    return pf.open(outname)

def divide(A,B, outname):
    A,B = gunzip(A,B)
    imarith(A, "/", B, outname)
    gzip(A,B)

    return pf.open(outname)


def combines(A,B,C,D, outname):
    ''' Creates outname with A+B+C+D '''
    if os.path.exists(outname):
        return pf.open(outname)

    try: os.remove('AB.fits')
    except: pass
    try: os.remove('CD.fits')
    except: pass
    add(A,B, 'AB.fits')
    add(C,D, 'CD.fits')

    add('AB.fits', 'CD.fits', outname)
    try: os.remove('AB.fits')
    except: pass
    try: os.remove('CD.fits')
    except: pass


def combine4(A,B,C,D, outname):
    '''Creates outname which is == A-(B+C+D)/3'''

    if os.path.exists(outname):
        return pf.open(outname)

    try: os.remove("APB.fits")
    except: pass
    add(A,B,"APB.fits")
    try: os.remove("3ABC.fits")
    except: pass
    add(B,C,"3ABC.fits")
    try: os.remove(outname)
    except: pass
    divide("3ABC.fits", "3", "ABC.fits")
    try: os.remove("3ABC.fits")
    except: pass
    try: os.remove(outname)
    except: pass
    subtract(A, "ABC.fits", outname)
    try: os.remove("ABC.fits")
    except: pass
    try: os.remove("APB.fits")
    except: pass

    return pf.open(outname)

def bgd_level(extractions):
    '''Remove background from extractions'''

    levels = []
    for spectrum in extractions:
        if spectrum.__dict__.has_key('spec') and spectrum.spec is not None \
            and spectrum.lamcoeff is not None:
        
            l, Fl = spectrum.get_counts(the_spec='specw')

            levels.append(np.median(Fl))

    bgd = np.median(levels)
    sd = np.std(levels)
    pl.plot(levels,'x')
    pl.axhline(bgd)
    pl.axhline(bgd+sd)
    pl.axhline(bgd-sd)
    pl.ylim(-20*sd-bgd,20*sd-bgd)
    pl.show()


def handle_extract(data, outname=None, fine='fine.npy',flexure_x_corr_nm=0.0,
    flexure_y_corr_pix = 0.0):

    exfile = "extracted_%s.npy" % outname
    if not os.path.exists(outname + ".npy"):
        E = Wavelength.wavelength_extract(data, fine, filename=outname,
            flexure_x_corr_nm = flexure_x_corr_nm,
            flexure_y_corr_pix= flexure_y_corr_pix,
            flat_corrections=flat_corrections)
            
        np.save(exfile, [E, meta])
    else:
        E, meta = np.load(exfile)

    return E

def handle_A(A, fine, outname=None, standard=None, corrfile=None,
    Aoffset=None, radius=2, flat_corrections=None):
    '''Loads 2k x 2k IFU frame "A" and extracts spectra from the locations
    in "fine". 

    Args:
        A (string): filename of ifu FITS file to extract from.
        fine (string): filename of NumPy file with locations + wavelength
            soln
        outname (string): filename to write results to
        Aoffset (2tuple): X (nm)/Y (pix) shift to apply for flexure correction
        radius (float): Extraction radius in arcsecond
        flat_corrections (list): A list of FlatCorrection objects for
            correcting the extraction

    Returns:
        The extracted spectrum, a dictionary:
        {'ph_10m_nm': Flux in photon / 10 m / nanometer integrated
        'nm': Wavelength solution in nm
        'N_spax': Total number of spaxels that created ph_10m_nm
        'skyph': Sky flux in photon / 10 m / nanometer / spaxel
        'radius_as': Extraction radius in arcsec
        'pos': X/Y extraction location of spectrum in arcsec}

    Raises:
        None
    '''

    fine = np.load(fine)
    if outname is None:
        outname = "%s" % (A)

    spec = pf.open(A)

    if Aoffset is not None:
        ff = np.load(Aoffset)
        flexure_x_corr_nm = ff[0]['dXnm']
        flexure_y_corr_pix = ff[0]['dYpix']
    else:
        flexure_x_corr_nm = 0
        flexure_y_corr_pix = 0

    if os.path.isfile(outname+".npy"):
        print "USING extractions in %s!" % outname
        print "rm %s.npy # if you want to recreate extractions" % outname
        E, meta = np.load(outname+".npy")
    else:
        print "CREATING extractions ..."
        E, meta = Wavelength.wavelength_extract(spec, fine, filename=outname,
            flexure_x_corr_nm=flexure_x_corr_nm, 
            flexure_y_corr_pix=flexure_y_corr_pix,
            flat_corrections = flat_corrections)

        np.save(outname, [E, meta])

    meta['airmass'] = spec[0].header['airmass']
    six, pos = identify_spectra_gui(E, radius=radius)
    skyix = identify_bgd_spectra(E, pos, inner=radius*1.1)
    res = interp_spectra(E, six, outname=outname+".pdf", corrfile=corrfile)
    sky = interp_spectra(E, skyix, onto=res[0]['nm'], outname=outname+"_sky.pdf", corrfile=corrfile)
    
    to_image(E, meta, outname, posA=pos)
    if standard is not None:
        print "STANDARD"
        wav = standard[:,0]/10.0
        flux = standard[:,1]

        fun = interp1d(wav, flux, bounds_error=False, fill_value = np.nan)
        correction = fun(res[0]['nm'])/res[0]['ph_10m_nm']

        res[0]['std-correction'] = correction


    airmass = meta['airmass']
    extCorr = 10**(Atm.ext(res[0]['nm']*10) * airmass/2.5)
    print "Median airmass corr: ", np.median(extCorr)

    ff = interp1d(sky[0]['nm'], sky[0]['ph_10m_nm'], bounds_error=False)
    skybgd = ff(res[0]['nm'])

    res[0]['exptime'] = spec[0].header['exptime']
    res[0]['Extinction Correction'] = 'Applied using Hayes & Latham'
    res[0]['extinction_corr'] = extCorr
    res[0]['skynm'] = sky[0]['nm']
    res[0]['skyph'] = sky[0]['ph_10m_nm']

    res[0]['ph_10m_nm'] -= skybgd 
    res[0]['ph_10m_nm'] *= extCorr * len(six)

    res[0]['radius_as'] = radius
    res[0]['position'] = pos
    res[0]['N_spax'] = len(six)
    res[0]['meta'] = meta
    res[0]['object_spaxel_ids'] = six
    res[0]['sky_spaxel_ids'] = skyix
    res[0]['sky_spectra'] = sky[0]['spectra']

    np.save("sp_" + outname, res)



def handle_AB(A, B, fine, outname=None, corrfile=None,
    Aoffset=None, Boffset=None, radius=2, flat_corrections=None):
    '''Loads 2k x 2k IFU frame "A" and "B" and extracts A-B and A+B spectra
    from the "fine" location. 

    Args:
        A (string): filename of ifu FITS file to extract from.
        B (string): filename of ifu FITS file to extract from.
        fine (string): filename of NumPy file with locations + wavelength
            soln
        outname (string): filename to write results to
        Aoffset (2tuple): X (nm)/Y (pix) shift to apply for flexure correction
        Boffset (2tuple): X (nm)/Y (pix) shift to apply for flexure correction
        radius (float): Extraction radius in arcsecond
        flat_corrections (list): A list of FlatCorrection objects for
            correcting the extraction

    Returns:
        The extracted spectrum, a dictionary:
        {'ph_10m_nm': Flux in photon / 10 m / nanometer integrated
        'var'
        'nm': Wavelength solution in nm
        'N_spaxA': Total number of "A" spaxels 
        'N_spaxB': Total number of "B" spaxels
        'skyph': Sky flux in photon / 10 m / nanometer / spaxel
        'radius_as': Extraction radius in arcsec
        'pos': X/Y extraction location of spectrum in arcsec}

    Raises:
        None
    '''

    fine = np.load(fine)
    if outname is None:
        outname = "%sm%s" % (A,B)

    if Aoffset is not None:
        ff = np.load(Aoffset)
        flexure_x_corr_nm = ff[0]['dXnm']
        flexure_y_corr_pix = ff[0]['dYpix']
    else:
        flexure_x_corr_nm = 0
        flexure_y_corr_pix = 0

    if os.path.isfile(outname + ".fits.npy"):
        print "USING extractions in %s!" % outname
        E, meta = np.load(outname + ".fits.npy")
        E_var, meta_var = np.load("var_" + outname + ".fits.npy")
    else:
        if not outname.endswith(".fits"): 
            outname = outname + ".fits"
            diff = subtract(A,B, outname)
            add(A,B, "tmpvar_" + outname)

            adcspeed = diff[0].header["ADCSPEED"]
            if adcspeed == 2: read_var = 22*22
            else: read_var = 5*5

        var = add("tmpvar_" + outname, str(read_var), "var_" + outname)
        os.remove("tmpvar_" + outname + ".gz")


        E, meta = Wavelength.wavelength_extract(diff, fine, 
            filename=outname,
            flexure_x_corr_nm = flexure_x_corr_nm, 
            flexure_y_corr_pix = flexure_y_corr_pix,
            flat_corrections=flat_corrections)
        meta['airmass1'] = diff[0].header['airmass1']
        meta['airmass2'] = diff[0].header['airmass2']
        meta['exptime'] = diff[0].header['exptime']
        np.save(outname, [E, meta])

        exfile = "extracted_var_%s.npy" % outname
        E_var, meta_var = Wavelength.wavelength_extract(var, fine, 
            filename=outname,
            flexure_x_corr_nm = flexure_x_corr_nm, 
            flexure_y_corr_pix = flexure_y_corr_pix,
            flat_corrections=flat_corrections)

        np.save("var_" + outname, [E_var, meta_var])

    sixA, posA = identify_spectra_gui(E, radius=radius)
    sixB, posB = identify_spectra_gui(E, radius=radius)

    to_image(E, meta, outname, posA=posA, posB=posB)

    skyA = identify_bgd_spectra(E, posA)
    skyB = identify_bgd_spectra(E, posB)

    allix = np.concatenate([sixA, sixB])
    resA = interp_spectra(E, sixA, sign=1, outname=outname+"_A.pdf", corrfile=corrfile)
    resB = interp_spectra(E, sixB, sign=-1, outname=outname+"_B.pdf", corrfile=corrfile)
    skyA = interp_spectra(E, skyA, sign=1, outname=outname+"_skyA.pdf", corrfile=corrfile)
    skyB = interp_spectra(E, skyB, sign=-1, outname=outname+"_skYB.pdf", corrfile=corrfile)
    varA = interp_spectra(E_var, sixA, sign=1, outname=outname+"_A_var.pdf", corrfile=corrfile)
    varB = interp_spectra(E_var, sixB, sign=1, outname=outname+"_B_var.pdf", corrfile=corrfile)
    
    ## Plot out the X/Y selected spectra
    XSA = []
    YSA = []
    XSB = []
    YSB = []
    for ix in sixA:
        XSA.append(E[ix].X_as)
        YSA.append(E[ix].Y_as)
    for ix in sixB:
        XSB.append(E[ix].X_as)
        YSB.append(E[ix].Y_as)

    pl.figure()
    pl.clf()
    pl.ylim(-30,30)
    pl.xlim(-30,30)
    pl.scatter(XSA,YSA, color='blue', marker='H', linewidth=.1)
    pl.scatter(XSB,YSB, color='red', marker='H', linewidth=.1)
    pl.savefig("XYs_%s.pdf" % outname)
    pl.close()
    # / End Plot

    np.save("sp_A_" + outname, resA)
    np.save("sp_B_" + outname, resB)
    np.save("var_A_" + outname, varA)
    np.save("var_B_" + outname, varB)

    ll = Wavelength.fiducial_spectrum()
    sky_A = interp1d(skyA[0]['nm'], skyA[0]['ph_10m_nm'], bounds_error=False)
    sky_B = interp1d(skyB[0]['nm'], skyB[0]['ph_10m_nm'], bounds_error=False)

    sky = np.nanmean([sky_A(ll), sky_B(ll)], axis=0)

    res = np.copy(resA)
    res = [{"doc": resA[0]["doc"], "ph_10m_nm": np.copy(resA[0]["ph_10m_nm"]),
        "nm": np.copy(resA[0]["ph_10m_nm"])}]
    res[0]['nm'] = np.copy(ll)
    f1 = interp1d(resA[0]['nm'], resA[0]['ph_10m_nm'], bounds_error=False)
    f2 = interp1d(resB[0]['nm'], resB[0]['ph_10m_nm'], bounds_error=False)

    airmassA = meta['airmass1']
    airmassB = meta['airmass2']

    extCorrA = 10**(Atm.ext(ll*10)*airmassA/2.5)
    extCorrB = 10**(Atm.ext(ll*10)*airmassB/2.5)
    print "Median airmass corr: ", np.median(extCorrA), np.median(extCorrB)
    res[0]['ph_10m_nm'] = \
        np.nansum([
            (f1(ll)-sky_A(ll)) * extCorrA, 
            (f2(ll)-sky_B(ll)) * extCorrB], axis=0) * \
            (len(sixA) + len(sixB))

    res[0]['exptime'] = meta['exptime']
    res[0]['Extinction Correction'] = 'Applied using Hayes & Latham'
    res[0]['extinction_corr_A'] = extCorrA
    res[0]['extinction_corr_B'] = extCorrB
    res[0]['skyph'] = sky
    res[0]['var'] = np.nanmean([varA[0]['ph_10m_nm'], varB[0]['ph_10m_nm']], 0) * (len(sixA) + len(sixB))
    res[0]['radius_as'] = radius
    res[0]['positionA'] = posA
    res[0]['positionB'] = posA
    res[0]['N_spaxA'] = len(sixA)
    res[0]['N_spaxB'] = len(sixB)
    res[0]['meta'] = meta
    res[0]['object_spaxel_ids_A'] = sixA
    res[0]['sky_spaxel_ids_A'] = skyA 
    res[0]['object_spaxel_ids_B'] = sixB
    res[0]['sky_spaxel_ids_B'] = skyB

    coef = chebfit(np.arange(len(ll)), ll, 4)
    xs = np.arange(len(ll)+1)
    newll = chebval(xs, coef)

    res[0]['dlam'] = np.diff(newll)

    np.save("sp_" + outname, res)

def measure_flexure(sky):
    ''' Measure expected (589.3 nm) - measured emission line in nm'''
    ll, ss = sky['nm'], sky['ph_10m_nm']

    pl.figure()
    pl.step(ll, ss)

    pix = SG.argrelmax(ss, order=20)[0]
    skynms = chebval(pix, sky['coefficients'])
    for s in skynms: pl.axvline(s)

    ixmin = np.argmin(np.abs(skynms - 589.3))
    dnm = 589.3 - skynms[ixmin]

    return dnm


def handle_ABCD(A, B, C, D, fine, outname=None, nsiglo=-1, nsighi=1, corrfile=None,flexure_pix_corr=(0,0)):
    '''Processes A-B files.

    1. Subtracts files to create A-(B+C+D), B-(A+C+D), C-(A+B+D), D-(A+B+C)
    2. Extracts the spectra from the four files
    3. Identifies the object as some multiple of the sigma
    
    Args:
        [ABCD]: string of filename
        fine: data structure containing the fine wavelength solution
        outname: Output file name, default is not to write output
        nsig[lo|hi]: float of the number of sigma to extract from (lo|hi)
        corrfile: Atmospheric extinction correction filename
        flexure_correction: 2-tuple containing float of pixel shift to apply in X, Y to correct flexure

    '''


    raise Exception("This code is not ready for prime time")

    if offset is None: offset = 0
    fine = np.load(fine)
    if outname is None: raise Exception("need output name")

    Ad = combine4(A, B,C,D, '%s_BCD.fits' % outname) 
    Bd = combine4(B, A,C,D, '%s_ACD.fits' % outname) 
    Cd = combine4(C, A,B,D, '%s_ABD.fits' % outname) 
    Dd = combine4(D, A,B,C, '%s_ABC.fits' % outname) 
    Sky = combines(A,B,C,D, '%s_ABCD.fits' % outname)

    Ae = handle_extract(Ad, outname='%s_BCD.fits' % outname, fine=fine, 
        offset=offset) 
    Be = handle_extract(Bd, outname='%s_ACD.fits' % outname, fine=fine,
        offset=offset) 
    Ce = handle_extract(Cd, outname='%s_ABD.fits' % outname, fine=fine,
        offset=offset) 
    De = handle_extract(Dd, outname='%s_ABC.fits' % outname, fine=fine,
        offset=offset) 
    Se = handle_extract(Sky, outname='%s_ABCD.fits' % outname, fine=fine,
        offset=offset) 


    sixA, posA = identify_spectra_gui(Ae, outname=outname+"_A.pdf")
    sixB, posB = identify_spectra_gui(Be, outname=outname+"_B.pdf")
    sixC, posC = identify_spectra_gui(Ce, outname=outname+"_C.pdf")
    sixD, posD = identify_spectra_gui(De, outname=outname+"_D.pdf")

    six1 = np.array(sixA.tolist() + sixB.tolist() )
    six2 = np.array(sixC.tolist() + sixD.tolist())

    sky1 = interp_spectra(Se, six1, outname=outname+"_sky.pdf", 
        corrfile=corrfile)[0]

    dnm = 0
    sky2 = interp_spectra(Se, six2, outname=outname+"_sky.pdf", 
        corrfile=corrfile, dnm=dnm)[0]

    print "CHECK: corrected sky spectra are offset by %s nm" % check

    resA = interp_spectra(Ae, sixA, outname=outname+"_A.pdf", 
        corrfile=corrfile, dnm=dnm)[0]
    resB = interp_spectra(Be, sixB, outname=outname+"_B.pdf", 
        corrfile=corrfile, dnm=dnm)[0]
    resC = interp_spectra(Ce, sixC, outname=outname+"_C.pdf", 
        corrfile=corrfile, dnm=dnm)[0]
    resD = interp_spectra(De, sixD, outname=outname+"_D.pdf", 
        corrfile=corrfile,dnm=dnm)[0]

    np.save(outname + "A.npy", [resA])
    np.save(outname + "B.npy", [resB])
    np.save(outname + "C.npy", [resC])
    np.save(outname + "D.npy", [resD])
    l = resA['nm']
    sp = resA['ph_10m_nm']
    cc  = resA['corrected-spec']

    fB = interp1d(resB['nm'], resB['ph_10m_nm'], fill_value=0,
        bounds_error=False)
    fC = interp1d(resC['nm'], resC['ph_10m_nm'], fill_value=0, 
        bounds_error=False)
    fD = interp1d(resD['nm'], resD['ph_10m_nm'], fill_value=0,
        bounds_error=False)

    cB = interp1d(resB['nm'], resB['corrected-spec'], fill_value=0,
        bounds_error=False)
    cC = interp1d(resC['nm'], resC['corrected-spec'], fill_value=0, 
        bounds_error=False)
    cD = interp1d(resD['nm'], resD['corrected-spec'], fill_value=0,
        bounds_error=False)

    sp += fB(l) + fC(l) + fD(l)
    cc += cB(l) + cC(l) + cD(l)


    result = [{"nm": l, "ph_10m_nm": sp, "corrected-spec": cc}]

    np.save(outname + ".npy", result)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=\
        '''Extracter.py:

            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--A', type=str, help='FITS A file')
    parser.add_argument('--B', type=str, help='FITS B file')
    parser.add_argument('--C', type=str, help='FITS C file')
    parser.add_argument('--D', type=str, help='FITS B file')
    parser.add_argument('fine', type=str, help='Numpy fine wavelength solution')
    parser.add_argument('--outname', type=str, help='Prefix output name')
    parser.add_argument('--nsiglo', type=float, help='Number sigma to extract below', default=-2)
    parser.add_argument('--nsighi', type=float, help='Number sigma to extract above', default=2)
    parser.add_argument('--std', type=str, help='Name of standard')
    parser.add_argument('--correction', type=str, help='Name of atmospheric correction file')
    parser.add_argument('--Aoffset', type=str, help='Name of "A" file that holds flexure offset correction information')
    parser.add_argument('--Boffset', type=str, help='Name of "B" file that holds flexure offset correction information')
    parser.add_argument('--Coffset', type=str, help='Name of "C" file that holds flexure offset correction information')
    parser.add_argument('--Doffset', type=str, help='Name of "D" file that holds flexure offset correction information')
    parser.add_argument('--radius_as', type=float, help='Extraction radius in arcsecond', default=3)
    parser.add_argument('--flat_correction', type=str, help='Name of flat field .npy file', default=None)

    args = parser.parse_args()


    if args.outname is not None:
        args.outname = args.outname.rstrip('.npy')

    if args.flat_correction is not None:
        print "Using flat data in %s" % args.flat_correction
        flat = np.load(args.flat_correction)
    else: flat = None

    if args.A is not None and args.B is not None and args.C is not None \
        and args.D is not None:
        handle_ABCD(args.A, args.B, args.C, args.D, args.fine, 
                outname=args.outname,
            nsiglo=args.nsiglo, nsighi=args.nsighi, corrfile=args.correction,
            offset=args.Aoffset, flat_corrections=flat)

    elif args.A is not None and args.B is not None:
        print "Handle AB"
        handle_AB(args.A, args.B, args.fine, outname=args.outname,
            corrfile=args.correction,
            Aoffset=args.Aoffset, Boffset=args.Boffset, 
            radius=args.radius_as, flat_corrections=flat)

    elif args.A is not None:
        if args.std is None:
            handle_A(args.A, args.fine, outname=args.outname,
                corrfile=args.correction,
                Aoffset=args.Aoffset, radius=args.radius_as,
                flat_corrections=flat)
        else:
            star = Stds.Standards[args.std]
            handle_A(args.A, args.fine, outname=args.outname,
                standard=star,
                Aoffset=args.Aoffset, radius=args.radius_as,
                flat_corrections=flat)
            
    else:
        print "I do not understand your intent, you must specify --A, at least"
