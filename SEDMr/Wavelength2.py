import argparse
import pdb
from multiprocessing import Pool
import numpy as np
import pylab as pl
import pyfits as pf
import scipy
import sys


import NPK.Fit 
import NPK.Bar as Bar
from astropy.table import Table 

from scipy.spatial import KDTree 
import scipy.signal as SG

from numpy.polynomial.chebyshev import chebfit, chebval

import SEDMr.Extraction as Extraction
from scipy.interpolate import interp1d
reload(NPK.Fit)
reload(Extraction)

from numpy import NaN, Inf, arange, isscalar, asarray, array
 
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
 
    return array(maxtab), array(mintab)



def read_catalog(catname):
    '''Read a sextractor catalog called catname and return'''


    cat = Table.read(catname, format='ascii.sextractor')

    return cat

def hg_to_kdtree(assoc_hg_spec):
    '''Convert mercury association to a KD Tree'''

    xs = []
    ys = []
    ids = []
    for id,f in enumerate(assoc_hg_spec):
        if not f.has_key(546.1): continue
        xs.append(f[546.1])
        ff = np.poly1d(f['coeff_ys'])
        ys.append(ff(xs[-1]))
        ids.append(f['seg_cnt'])

    
    data = np.array([xs,ys])
    return KDTree(data.T), np.array(ids)

def fiducial_spectrum(lamstart=1000.0, lamratio=239./240., len=250):
    '''Return a typical SED Machine spectrum, use for interpolating on grid

                                        x
    Equation looks like 1000 x (239/240)

    Args:
        lamstart(float): default is 1000 nm
        lamratio(float): default is 239/240, good for SEDM
        len(int): default of 250 yields a spectrum that goes to ~ 350 nm'''

    xx = np.arange(len)

    return lamstart * lamratio**xx

def assoc_hg_with_flats(domedat, hgcat, guess_offset= {365.0: 231,
    404.6:214, 435.8:193, 546.1: 133, 578.00: 110}, outname='assoc_Hg'):
    
    '''Given a set of functions defining the ridgeline of dome flats
    and a list of positions of mercury lamps, provide a crude wavelength
    solution

    Args:
        guess_offset: Dictionary of {wavelength in nm: x position} indicating
            the rough distance from the segment start to the wavelength
    '''

    wavetrees = {}
    for k,v in hgcat.iteritems():
        wavetrees[k] = KDTree(v)

    reg = ''
    width = Bar.setup()
    update_rate = len(domedat)/width
    for idx,spec_pos in enumerate(domedat):
        if idx % update_rate == 0: Bar.update()
        if not spec_pos['ok']: continue
        tracefun = np.poly1d(spec_pos['coeff_ys'])
        minx = spec_pos['xs'][0]
        

        for wavelen,v in wavetrees.iteritems():
            offset = guess_offset[wavelen]
            pt = (minx + offset, tracefun(minx+offset))
            results = v.query_ball_point(pt, 15)

            for res in results:
                x_hg,y_hg = v.data[res]
                y_trace = tracefun(x_hg)

                if np.abs(y_hg - y_trace) < 3:
                    spec_pos[wavelen] = x_hg
                    if wavelen == 578.0:
                        reg += 'point(%s,%s) # point=x 5 text={%s}\n' % (x_hg, y_trace, spec_pos['seg_cnt'])
                    else:
                        reg += 'point(%s,%s) # point=x 5\n' % (x_hg, y_trace)

    Bar.done()
    f = open("%s.reg" % outname, "w")
    f.write(reg)
    f.close()

    np.save("%s.npy" % outname, [domedat])

    return spec_pos

        
def find_hg_spectra(lines, dYlimit=2, outname="find_spectra"):
    '''Based on line ratios and positions, determine positions of Hg lines
    
    Args:
        lines: line list
        dYlimit: 

    Results:
        Writes [return value] to outname in npy format
    
    Return:
        Returns dictionary of {wavelength nm: [(pixel coords) ...], ...}
        Each wavelength has a list of pixel coords of different length.
        The coordinates are diassociated from one and other.
        '''


    data = []

    for line in lines:
        data.append((line['X_IMAGE'], line['Y_IMAGE']))
    
    data = np.array(data)

    kdt = KDTree(data)

    reg = ""
    ratios = []
    hg546 = []
    hg578 = []
    hg436 = []
    hg405 = []
    hg365 = []
    for line in lines:
        point = [line['X_IMAGE'], line['Y_IMAGE']]
        results = kdt.query_ball_point(point, 38)
        X, Y = point
        num = line['NUMBER']
        flux = float(line['FLUX_ISO'])
        mean_flux = np.median(flux)/2

        for resix in results:
            resX, resY = data[resix]
            resFlux = float(lines[resix]['FLUX_ISO'])
            dist = np.abs(resX-X)
            dX = resX - X
            
            dY = np.abs(resY - Y)
            if dY > dYlimit: continue

            if resFlux < 500: continue

            ratio = flux/resFlux
            bright = (flux > mean_flux)
            resbright = (resFlux > mean_flux)

            if (2 < ratio < 7) and (-16 < dX < -12) and bright:
                reg += 'circle(%s,%s, 3) # color=red\n' % (X, Y)
                reg += 'circle(%s,%s, 3) # color=magenta\n' % (resX, resY)
                hg546.append((X,Y))
                hg578.append((resX,resY))
                continue

            if (1/6. < ratio < 2) and (-24 < dX < -18) and resbright:
                hg405.append((X,Y))
                hg436.append((resX,resY))
                reg += 'circle(%s,%s, 3) # color=green\n' % (X, Y)
                reg += 'circle(%s,%s, 3) # color=yellow\n' % (resX, resY)
                continue

            if (4 < ratio < 70) and (20 < dX < 38) and bright and (not resbright):
                hg365.append((resX,resY))
                reg += 'circle(%s,%s, 3) # color=blue\n' % (resX, resY)
                continue

    f = open(outname + ".reg", "w")
    f.write(reg)
    f.close()
    
    res = [{365.0: np.array(hg365),
            404.6: np.array(hg405),
            435.8: np.array(hg436),
            546.1: np.array(hg546),
            578.00: np.array(hg578)}]


    np.save(outname, res)

    return res[0]


def wavelength_extract(HDUlist, wavecalib, filename='extracted_spectra.npy',
    flexure_x_corr_nm = 0.0, flexure_y_corr_pix = 0.0, extract_width=3):
    
    dat = HDUlist[0].data
    exptime = HDUlist[0].header['EXPTIME']

    extractions = []
    update_rate = len(wavecalib) / Bar.setup()


    print "Applying %f/%f offset" % (flexure_x_corr_nm, flexure_y_corr_pix)
    flexure_y_corr_pix = np.round(flexure_y_corr_pix)
    for ix, ss in enumerate(wavecalib):
        if ix % update_rate == 0: Bar.update()
        if not ss.ok: 
            extractions.append(Extraction.Extraction(seg_id=ss.seg_id,
                ok=False, trace_sigma=ss.trace_sigma))
            continue
        minx = np.max((0,ss.xrange[0]-5))
        maxx = np.min((minx + 265,2047))
        yfun = np.poly1d(ss.poly)

        xpos = xrange(minx, maxx)
        res = np.zeros(len(xpos))
        res[:] = np.nan
        resw = np.zeros(len(xpos))
        resw[:] = np.nan
        sigma2 = ss.trace_sigma * ss.trace_sigma
        if np.isnan(sigma2): sigma2 = 4.

        for i in xrange(len(xpos)):
            X = xpos[i]
            Y = yfun(X)
            if not np.isfinite(X) or not np.isfinite(Y):
                continue

            # Extract width requires asymmetry in the slice
            # slice[-2:3] will return elements -2 to +2 around 0
            # e.g. len(slice[-2:3]) == 5
            Ys = slice(
                np.max((0,np.int(Y)+flexure_y_corr_pix-extract_width)),
                np.min((np.int(Y)+flexure_y_corr_pix+extract_width+1, 2047)))

            profile = np.arange(np.round(Ys.stop)-np.round(Ys.start))

            NNN = dat[Ys,X].shape[0]
            profile = profile[0:NNN]
            profile -= (len(profile)-1)/2.0
            profile = np.exp(- profile*profile/(2*sigma2))
            profile /= np.mean(profile)



            res[i] = np.sum(dat[Ys,X])
            resw[i]= np.sum(dat[Ys,X]*profile)


        extractions.append( 
            Extraction.Extraction(xrange=(minx,maxx), yrange=(yfun(xpos[0]), 
                                    yfun(xpos[-1])),
                                    poly=ss.poly, spec=res/exptime, 
                                    specw=resw/exptime, 
                                    seg_id=ss.seg_id, exptime=exptime, ok=True, 
                                    trace_sigma=ss.trace_sigma))

        if ss.__dict__.has_key('lamcoeff') and ss.lamcoeff is not None:
            extractions[-1].lamcoeff = ss.lamcoeff.copy()
            extractions[-1].lamcoeff[0] -= flexure_x_corr_nm
            extractions[-1].lamrms = ss.lamrms

        if ss.__dict__.has_key('mdn_coeff') and ss.mdn_coeff is not None:
            extractions[-1].mdn_coeff = ss.mdn_coeff.copy()
            extractions[-1].mdn_coeff[0] -= flexure_x_corr_nm
            extractions[-1].lamrms = ss.lamrms

    Bar.done()

    np.save(filename, extractions)
    return extractions

def extract_helper(ss):
    global dat

    if not ss['ok']: 
        return Extraction.Extraction(seg_id=ss['seg_cnt'], 
            trace_sigma=ss['trace_sigma'], ok=False)

    minx = np.max((0,ss['xs'][0]-5))
    maxx = np.min((minx + 265,2047))
    yfun = np.poly1d(ss['coeff_ys'])

    xpos = xrange(minx, maxx)
    res = np.zeros(len(xpos))
    res[:] = np.nan
    resw = np.zeros(len(xpos))
    resw[:] = np.nan
    sigma2 = ss['trace_sigma']*ss['trace_sigma']
    if np.isnan(sigma2): sigma2 = 4.

    for i in xrange(len(xpos)):
        X = xpos[i]
        Y = yfun(X)
        if Y < 0: Y = 0
        if Y > 2046: Y = 2046
        if not np.isfinite(X) or not np.isfinite(Y):
            continue
        Ys = slice(np.max((0,np.int(Y)-3)),np.min((np.int(Y)+3, 2047)))

        profile = np.arange(Ys.stop-Ys.start)
        profile -= (len(profile)-1)/2.0
        profile = np.exp(- profile*profile/(2*sigma2))
        profile /= np.mean(profile)
        res[i] = np.sum(dat[Ys,X])
        try:
            resw[i]= np.sum(dat[Ys,X]*profile)
        except:
            import pdb
            pdb.set_trace()

    hg_lines = {}
    for wave, pix in ss.iteritems():
        if type(wave) is float:
            hg_lines[wave] = pix-minx-1
    

    return Extraction.Extraction(xrange=(minx,maxx), yrange=(yfun(xpos[0]), 
                                yfun(xpos[-1])),
                                poly=ss['coeff_ys'], spec=res, specw=resw,
                                hg_lines=hg_lines, seg_id=ss['seg_cnt'],
                                ok=True, trace_sigma=ss['trace_sigma'])

    




def extract(HDUlist, assoc_hg_spec, filename='raw_extractions'):
    global dat


    dat = HDUlist[0].data

    
    extractions = []

    p = Pool(16)
    extractions = p.map(extract_helper, assoc_hg_spec)
    p.close()

    np.save(filename, extractions)
    return extractions


def median_fine_grid(fine):
    ''' Using the 7 nearest neighbors, median the wavelength solution. Refit the wavelength solution to this median wavelength.

    Input:
        fine: Dictionary that contains the fine wavelength solution

    Returns:
        fine.mdn_coeff contains updated chebyshev polynomial coefficients '''
        
    xs = []
    ys = []
    ids = []
    for gix, g in enumerate(fine):
        if not g.ok: 
            xs.append(-999)
            ys.append(-999)
            ids.append(None)
            continue

        xs.append(np.mean(g.xrange))
        ys.append(np.mean(g.yrange))
        ids.append(g.seg_id)

    xs,ys,ids = map(np.array, [xs,ys,ids])
    dat = np.array((xs,ys)).T
    dat[dat != dat] = -999 # Correct NaNs
    KD = KDTree(dat)

    assert(len(ids) == len(fine))
    ixs = np.arange(0, 265, 1)
    for idx, spec in enumerate(fine):
        seg_id = spec.seg_id
        loc = dat[seg_id-1]

        if ids[idx] is None: continue
        dist, nearest_ixs = KD.query(loc, k=10)

        lls = []
        num_in = 0
        shifts = []
        for nearest in nearest_ixs[1:]:
            if nearest is None: continue
            if fine[nearest].hgcoef is None: continue
            if fine[nearest].lamcoeff is None: continue
            if fine[nearest].lamrms > 0.2: continue
            if np.abs(fine[nearest].xrange[1] - fine[nearest].xrange[0]) < 50: continue

            xx = np.arange(265)
            ll = chebval(xx+fine[nearest].xrange[0], fine[nearest].lamcoeff)
            if np.any(np.diff(ll) > 0): continue

            lls.append(ll)
            num_in += 1
            if num_in > 30: break

        lls = np.array(lls, dtype=np.float)
        if len(lls) == 0:
            spec.mdn_coeff = spec.lamcoeff
            continue

        try: new_lls = scipy.stats.nanmedian(lls, 0)
        except:
            import pdb
            pdb.set_trace()

        diff = (lls - new_lls) / new_lls
        stds = scipy.stats.nanstd(diff, 0)
        bad = np.abs(diff/stds) > 4
        if bad.any():
            lls[bad] = np.nan
            new_lls = scipy.stats.nanmedian(lls, 0)

        diff = (lls - new_lls) / new_lls
        stds = scipy.stats.nanstd(diff, 0)
        bad = np.abs(diff/stds) > 2
        if bad.any():
            lls[bad] = np.nan
            new_lls = scipy.stats.nanmedian(lls, 0)

        new_lls = scipy.stats.nanmedian(lls, 0)

        spec.mdn_coeff = chebfit(np.arange(len(new_lls))+spec.xrange[0], 
            new_lls, 3)

    pl.figure(3)
    pl.clf()
    pl.xlim(360, 1100)
    pl.ioff()
    pl.figure(4)
    pl.clf()
    pl.xlim(360, 1100)
    pl.ioff()
    for g in fine:
        if g.mdn_coeff is None: continue
        if g.specw is None: continue
        if g.hgcoef is None: continue
        if len(g.specw) < 30: continue


        if g.xrange[0] < 100: continue
        if g.xrange[0] > 1900: continue
        if g.yrange[0] < 100: continue
        if g.yrange[0] > 1900: continue

        ix = np.arange(len(g.specw))
        pl.figure(3)
        ll = chebval(ix+g.xrange[0], g.mdn_coeff)
        pl.plot(ll, g.specw, '.')
        pl.figure(4)
        try: 
            ll = chebval(ix+g.xrange[0], g.lamcoeff)
            pl.plot(ll, g.specw, '.')
        except: pass
        
    pl.ion()
    pl.show()

    return fine


def median_rough_grid(gridded, Hg_E, outname='median_rough_wavelength.npy'):
    ''' Using the 7 nearest neighbors, median the wavelength solution coefficients '''

    this_code_is_of_very_little_use()

    '''
        median_rough grid was used to test the idea
        as of 20 jan 2015, the algorithm tested by this code does not
        work.
    '''
    xs = []
    ys = []
    ids = []
    for gix, g in enumerate(gridded):
        if not g.ok: continue
        if 546.1 not in g.hg_lines: continue

        xs.append(g.hg_lines[546.1] + g.xrange[0])
        ys.append(np.mean(g.yrange))
        ids.append(g.seg_id)

    xs,ys,ids = map(np.array, [xs,ys,ids])
    dat = np.array((xs,ys)).T
    KD = KDTree(dat)

    for ix, id in enumerate(ids):
        loc = dat[ix]
        dist, nearest_ixs = KD.query(loc, k=14)
        #print dist, nearest_ixs

        lls = []
        num_in = 0
        for nix, nearest in enumerate(ids[nearest_ixs]):
            #print nearest, dist[nix]
            #if dist[nix] > 100: continue
            if gridded[nearest-1].hgcoef is None: continue
            xx = np.arange(265)
            ll = chebval(xx, gridded[nearest-1].hgcoef)
            if np.any(np.diff(ll) > 0): continue

            lls.append(ll)
            num_in += 1
            if num_in > 30: break

        if len(lls) == 0:
            import pdb
            pdb.set_trace()
        lls = np.array(lls)
        if idx > 500:
            import pdb
            pdb.set_trace()
        try: 
            new_lls = scipy.stats.nanmedian(lls, 0)
            new_lls_std = scipy.stats.nanstd(lls, 0)
        except:
            import pdb
            pdb.set_trace()

        diff = (lls - new_lls) / new_lls
        bad = np.abs(diff) > .05
        lls[bad] = np.nan
        new_lls = scipy.stats.nanmedian(lls, 0)

        diff = (lls - new_lls) / new_lls
        bad = np.abs(diff) > .05
        lls[bad] = np.nan
        new_lls = scipy.stats.nanmedian(lls, 0)

        gridded[ix].mdn_hgcoef = chebfit(np.arange(len(new_lls)), new_lls, 4)
        import pdb
        pdb.set_trace()

    pl.figure(3)
    pl.clf()
    pl.xlim(360, 1100)
    pl.ioff()
    for g in gridded:
        if "mdn_hgcoef" not in g.__dict__: continue
        if g.specw is None: continue
        if g.hgcoef is None: continue
        if len(g.specw) < 30: continue

        ix = np.arange(len(g.specw))
        ll = chebval(ix, g.mdn_hgcoef)
        pl.plot(ll, g.specw, '.')
        
    pl.ion()
    pl.show()
    import pdb
    pdb.set_trace()

def scale_on_547(spec):
    if spec is None: return None
    if spec.spec is None: return None
    if spec.hgcoef is None: return None

    ix = np.arange(0, len(spec.specw), .1)
    ix_547 = np.argmin(np.abs(chebval(ix, spec.hgcoef) - 547.0))

    scale = np.arange(len(spec.specw)) - ix[ix_547]

    return scale


def stretch_fit_helper(specno):
    ''' Helper for Multiprocessing Pool'''
    global pix_coeffs, fid_ix, fs1, fs2, slopes, squares, specs, funs

    spec = specs[specno]
    if spec is None: return None

    xcs1 = np.zeros((len(slopes), len(squares)))
    xcs2 = np.zeros((len(slopes), len(squares)))

    s1f, s2f = funs[specno]
    for i, slope in enumerate(slopes):
        for j, square in enumerate(squares):
            # Step 4
            ns1 = s1f(fid_ix * slope + fid_ix**2 * square)
            ns2 = s2f(fid_ix * slope + fid_ix**2 * square)
            xcs1[i,j] = np.nansum(fs1*fs1*ns1*ns1)
            xcs2[i,j] = np.nansum(fs2*fs2*ns2*ns2)

    slope1, square1 = np.unravel_index(np.argmax(xcs1), xcs1.shape)
    slope2, square2 = np.unravel_index(np.argmax(xcs2), xcs2.shape)

    ix1 = fid_ix * slopes[slope1] + fid_ix**2 * squares[square1]
    ix2 = fid_ix * slopes[slope2] + fid_ix**2 * squares[square2]

    tofit = np.zeros(len(ix1))
    tofit[0:150] = ix2[0:150]
    tofit[150:] = ix1[150:]

    # Step 5
    fc = chebfit(fid_ix, tofit, 2)
    xc = (xcs1[slope1, square1] + xcs2[slope2, square2] ) / np.nansum(fs1)

    print "%4.0i = %1.3e : %1.3f %+1.2e, %1.3f %+1.2e" % (specno, xc,slopes[slope1], squares[square1], slopes[slope2], squares[square2])

    if False:
        pl.figure(4)
        pl.plot(fid_ix, s1f(tofit)+s2f(tofit))

        import pdb
        pdb.set_trace()


    return fc, xc


def get_stretched(fid_ix, coeffs, spec1, spec2=None):
    ''' Reinterpolate spec1 and spec2 onto the fiducial index set by stretching and scaling the spectra

    This is used after a call that looks something like:
        ll = np.load("Hg_ext.npy")
        lx = np.load("Xe_ext.npy")
        gg = rough_grid(ll)
        coeffs, fix = stretch_set(gg, lx)

        s1, s2 = get_stretched(fix, coeffs[5], ll[5], lx[5])
    
    Returns:
        s1, s2 [float(len(fid_ix))]: Spectra'''
    newix = chebval(fid_ix, coeffs)
    ix = scale_on_547(spec1)
    sf1 = interp1d(ix, spec1.specw, bounds_error=False, fill_value=0)
    if spec2 is None: sf2 = lambda x: 0.0
    else:
        sf2 = interp1d(ix, spec2.specw, bounds_error=False, fill_value=0)

    s1 = sf1(newix)
    s2 = sf2(newix)

    return s1/s1.mean()*spec1.specw.mean(), s2/s2.mean()*spec2.specw.mean()


def stretch_set(Hg_set, Xe_set, mscales=None):
    ''' Shift and stretch spectra so they have the same pixel index 
    
    
    Steps:
    1. Shift spectra without interpolation onto a common pixel grid with index 0 being a prominent Hg line. 
    2. Create a set of interpolating functions for each of the spectra in the set. These functions are called as flux(pixel).
    3. Brute force search for the best stretch and 2nd-order polynomial coefficient for each spectrum in the set, compared to the fiducial spectrum.
    4. Brute force will measure the above for the xenon and mercury spectra indepdentnly.
    5. Then we stitch together a master pixel shift based on a chebyshev fit of the xenon and mercury spectra. The chebyshev fit creates a function such that f(pixel) --> pixel. On this new grid, all spectra are close to the same.
    '''

    # Global is for interprocess comms
    global pix_coeffs, fid_ix, fs1, fs2, slopes, squares, specs, funs

    assert(len(Hg_set) == len(Xe_set))

    
    # Step 1
    gridded_set = []
    for index, hg in enumerate(Hg_set):
        gridded_set.append(None)
        xe = Xe_set[index]
        if hg.spec is None: continue
        if xe.spec is None: continue
        if hg.hgcoef is None: continue

        assert(len(hg.spec) == len(xe.spec))

        ix = np.arange(len(hg.specw))
        ix_547 = np.argmin(np.abs(chebval(ix, hg.hgcoef) - 547.0))
        if ix_547 > 180: continue
        if ix_547 < 100: continue

        scale = np.arange(len(hg.spec)) - ix_547

        gridded_set[-1] = (scale, hg.spec, xe.spec) 


    fiducial = gridded_set[len(gridded_set)/2]
    if fiducial is None or len(fiducial[0]) != 265: 
        fiducial = gridded_set[len(gridded_set)/2 + 1]
    if fiducial is None or len(fiducial[0]) != 265: 
        fiducial = gridded_set[len(gridded_set)/2 + 2]
    if fiducial is None or len(fiducial[0]) != 265: 
        fiducial = gridded_set[len(gridded_set)/2 + 3]
    if fiducial is None or len(fiducial[0]) != 265: 
        fiducial = gridded_set[len(gridded_set)/2 + 4]


    print len(gridded_set)

    # Step 2
    pl.figure(3)
    pl.clf()


    def interp_functions():
        funs = []
        for num, element in enumerate(gridded_set):
            funs.append(None)
            if element is None: continue

            ix, s1, s2 = element
            s1f = interp1d(ix, s1, fill_value=np.nan, bounds_error=False)
            s2f = interp1d(ix, s2, fill_value=np.nan, bounds_error=False)

            funs[-1] = (s1f, s2f)

        return funs
        
    
    funs = interp_functions()



    pix_coeffs = []
    fid_ix, fs1, fs2 = fiducial
    # Note 0.001 over 250 pixels is about a quarter pixel
    # max of 0.5e-3 over 150**2 is a handful of pixels.
    slopes = np.arange(0.95,1/.95,.001)
    squares = np.linspace(-1e-3, 1e-3, 15)

    specs = gridded_set
    p = Pool(16)
    results = p.map(stretch_fit_helper, range(len(specs)))
    p.close()
    
    pix_coeffs = []
    xcs = []
    for res in results:
        if res is None:
            pix_coeffs.append(None)
            xcs.append(None)
        else:
            pix_coeffs.append(res[0])
            xcs.append(res[1])


    return pix_coeffs, fiducial, xcs


# TODO -- REFIT Wavelengths given a good guess

linelist = {
    "He": [587.5, 667.8, 706.5],
    "Cd": [467.8, 479.9, 508.5],
    "Hg": [578, 546.1, 435.8, 404.6, 365],
    "Xe": [764, 828]
}

def snap_solution_into_place(PARS): 
    ''' Return new Chebyshev coefficients for best fit wavelength solution '''

    if PARS is None:
        return (None, None)
    ixs, coef, lamp_spec = PARS

    fitfun = NPK.Fit.gaussian5
    resfun = NPK.Fit.mpfit_residuals(fitfun)


    def fitlinelist(cc):
        lams = chebval(ixs, cc)
        

        #pl.figure(1)
        #pl.clf()
        results = []
        wresults = []
        for lampname, lampspec in lamp_spec.iteritems():
            #pl.step(ixs, lampspec, where='mid')

            for line in linelist[lampname]:
                
                lix = np.argmin(np.abs(line-lams))

                if lix < 5 or lix > (len(ixs)-5): continue

                cutout = slice(lix-4, lix+4)
                xct = ixs[cutout]
                sct = lampspec[cutout]

                parguess = [
                    {'value': sct.max()-sct.min()},  # Scale
                    {'value': xct[sct.argmax()]}, # Centroid!!
                    {'value': 1.2}, # Sigma
                    {'value': sct.min()},# offset
                    {'value': 0}]  # slope

                pars = NPK.Fit.mpfit_do(resfun, xct, sct, parguess,
                    error=np.sqrt(sct))
                if pars.status != 1: continue

                results.append((line, pars.params[1], pars.perror[1]))
            

        results = np.array(results)

        if len(results) < 3:
            return cc, np.nan

        LS = results[:,0]
        IXS = results[:,1]
        WS = results[:,2]

        #for IX in IXS: pl.axvline(IX,color='red')

        nc = len(LS)-1
        if nc > 5: nc = 4
        newcoef = chebfit(IXS, LS, nc, w=WS)
        res = np.abs(chebval(IXS, newcoef) - LS)/LS

        if res[LS==365] < 200: 
            sl = np.where(LS==365)[0]
            np.delete(IXS, sl)
            np.delete(LS, sl)

            newcoef = chebfit(IXS, LS, nc, w=WS)
            res = np.abs(chebval(IXS, newcoef) - LS)/LS

        return newcoef,res


    newcoef,res = fitlinelist(coef)
    newcoef,newres = fitlinelist(newcoef)
    newcoef,newres = fitlinelist(newcoef)
    #pl.figure(2)
    #pl.clf()

    #pl.figure(2)
    #pl.plot(chebval(ixs, newcoef), lamp_spec['Hg'])
    #for line in linelist["Hg"]: pl.axvline(line, color='red')

    return newcoef, np.sqrt(np.mean(res*res))

def snap_solution_into_place_all(fine, Hgs, Xes, Cds=None, Hes=None):
    
    PARS = []
    for i in xrange(len(fine)):
        PARS.append(None)
        
        if fine[i].xrange is None: continue
        if fine[i].mdn_coeff is None: continue

        ixs = np.arange(*fine[i].xrange)
        coef = fine[i].mdn_coeff
        lamp_spec = {"Hg": fine[i].specw, 
            "Xe": Xes[i].specw}

        if Cds is not None: lamp_spec["Cd"] = Cds[i].specw
        if Hes is not None: lamp_spec["He"] = Hes[i].specw

        PARS[-1] = (ixs, coef, lamp_spec) 


    p = Pool(16)
    results = p.map(snap_solution_into_place, PARS)
    p.close()

    for i, res in enumerate(results):
        fine[i].lamcoeff = res[0]
        fine[i].lamrms = res[1]
        

    return fine
        


def fit_he_lines(SS, guesses = {587.5: -18, 667.8: -48, 706.5: -60}, plot=False,
    Ncutout=7):
    return fit_known_lines(SS, guesses, plot, Ncutout)

def fit_cd_lines(SS, guesses = {467.8: 41, 479.9: 33, 508.5: 18}, plot=False, 
    Ncutout=5):
    return fit_known_lines(SS, guesses, plot, Ncutout)

def fit_hg_lines(SS, guesses = {578: -14, 546.1: 0, 435.8: 60, 404.6: 81}, 
    plot=False, Ncutout=7):
    return fit_known_lines(SS, guesses, plot, Ncutout)

def fit_xe_lines(SS, guesses = {764: -77, 828: -93}, plot=False, Ncutout=5):
    return fit_known_lines(SS, guesses, plot, Ncutout)

def fit_known_lines(SS, guesses, plot, Ncutout):
    ''' Fit lines based on guess pixel positions against a fiducial spectrum 

    This function is mapable
    
    Args:
        SS is a list of two elements containting:
            spix(int[Ns]): Index values (pixel) of spec
            spec(float[Ns]): Flux from corresponding index values
        guesses: Dictionary containing {wavelength in nm: pixel guess position}
    Returns:
        [ (wavelength in nm, centroid position in pixel) ] with length 
            len(guesses.keys()).

        Intent is for this list to be fit with polynomials to construct the
        wavelength solution.
        
    '''

    if SS is None: return None
    spix, spec = SS

    res = []
    fitfun = NPK.Fit.gaussian5
    resfun = NPK.Fit.mpfit_residuals(fitfun)


    for lam, guesspos in guesses.iteritems():
        
        cutout = np.where(np.abs(spix - guesspos) < Ncutout)[0]
        xct = spix[cutout]
        sct = spec[cutout]

        parguess = [
            {'value': sct.max()-sct.min()},  # Scale
            {'value': xct[sct.argmax()]}, # Centroid!!
            {'value': 1.3}, # Sigma
            {'value': sct.min()},# offset
            {'value': 0}]  # slope

        pars = NPK.Fit.mpfit_do(resfun, xct, sct, parguess)
        if pars.status != 1: continue
        res.append([lam, pars.params[1]])

        if plot:
            pl.plot(xct, sct, 'o')
            pl.plot(xct, fitfun(pars.params, xct))

    if plot:
        import pdb
        pdb.set_trace()
        pl.clf()
    res = np.array(res)
    return res



def fit_all_lines(fiducial, hg_spec, xe_spec, xxs, cd_spec=None, he_spec=None):
    '''Fit mercury + xenon lines to a set of spectra that are put on a common fiducial grid.

    Args:
        fiducial(int[Nf]): Fiducial index values range from about -130 to + 130
        hg_spec(Extraction[Ne]): Mercury spectra
        xe_spec(Extraction[Ne]): Xenon spectra
        xxs(Gridded[Ne]): Results from stretch_set() call

    '''

    assert(len(hg_spec) > 500)
    assert(len(hg_spec) == len(xe_spec))
    if cd_spec is not None:
        assert(len(hg_spec) == len(cd_spec))
    if he_spec is not None:
        assert(len(hg_spec) == len(he_spec))

    assert(len(xxs) == len(hg_spec))


    Hgs = []
    Xes = []
    Cds = []
    Hes = []
    print "Stretching spectra"
    for i in xrange(len(hg_spec)):
        if hg_spec[i] is None or xxs[i] is None:
            Hgs.append(None)
            Xes.append(None)
            Cds.append(None)
            Hes.append(None)
            continue

        s1, s2 = get_stretched(fiducial, xxs[i], hg_spec[i],
                                spec2=xe_spec[i])

        if cd_spec is not None:
            s1, s3 = get_stretched(fiducial, xxs[i], hg_spec[i],
                                spec2=cd_spec[i])
            Cds.append( (fiducial, s3) )
        
        if he_spec is not None:
            s1, s4 = get_stretched(fiducial, xxs[i], hg_spec[i],
                                spec2=he_spec[i])
            Hes.append( (fiducial, s4) )

        Hgs.append( (fiducial, s1) )
        Xes.append( (fiducial, s2) )


    p = Pool(16)
    print "Fitting Hg lines"
    hg_locs = p.map(fit_hg_lines, Hgs)
    p.close()
    p = Pool(16)
    print "Fitting Xe lines"
    xe_locs = p.map(fit_xe_lines, Xes)
    p.close()
    
    if cd_spec is not None:
        p = Pool(16)
        print "Fitting Cd lines"
        cd_locs = p.map(fit_cd_lines, Cds)
        p.close()

    if he_spec is not None:
        p = Pool(16)
        print "Fitting He lines"
        he_locs = p.map(fit_he_lines, Hes)
        p.close()

    print "Determining best fit to all lines"
    fits = []
    residuals = []
    rss = []
    pl.figure(3)
    pl.clf()
    for i in xrange(len(hg_spec)):
        fits.append(None)
        rss.append(None)
        hg = hg_locs[i]
        xe = xe_locs[i]

        if cd_spec is not None: cd = cd_locs[i]
        if he_spec is not None: he = he_locs[i]

        if hg is None: continue
        if xe is None: continue
        if cd_spec is not None and cd is None: continue
        if he_spec is not None and he is None: continue

        try:
            ll = np.concatenate([hg[:,0], xe[:,0]])
            ix = np.concatenate([hg[:,1], xe[:,1]])

            if cd_spec is not None:
                ll = np.concatenate([ll, cd[:,0]])
                ix = np.concatenate([ix, cd[:,1]])
            if he_spec is not None:
                ll = np.concatenate([ll, he[:,0]])
                ix = np.concatenate([ix, he[:,1]])

        except:
            continue
        weights = np.ones(len(ll))
        weights[ll<400] = 0.03
        weights[ll>860] = 0.1
        

        if i == 50: 
            print ix
            print ll

        nc=3
        if len(ll) < 5: nc=2
        ff = chebfit(ix, ll, nc, w=weights)

        fits[-1] = ff

        res = {}
        for jx, l in enumerate(ll): res[l] = chebval(ix[jx], ff) - l
        residuals.append(res)

        rss[-1] = np.sqrt(np.sum((chebval(ix, ff) - ll)**2))

        if False:
            pl.figure(3)
            pl.step(chebval(fiducial, ff), Hgs[i][1]+Hes[i][1]+Xes[i][1])
            pl.figure(3)
            pl.step(chebval(fiducial, ff), Hgs[i][1]+Hes[i][1]+Xes[i][1])


    pl.clf()
    XS = []
    YS = []
    CS = []
    for i in xrange(len(hg_spec)):
        try:
            x = np.mean(hg_spec[i].xrange)
            y = np.mean(hg_spec[i].yrange)
            c  = rss[i]
        except:
            XS.append(None)
            YS.append(None)
            continue

        XS.append(x)
        YS.append(y)
        CS.append(c)


    return np.array(fits), residuals, np.array(rss), (XS, YS)

def coeffs_to_spec(fix, gridded, rgrd_coef, lam_coef):
    ''' Returns the spectrum given the unadultered spectrum and fits.

    Args:
        fix(array) -- Fiducial index positions
        gridded(Extraction) -- Extracted spectrum to plot
        rgrd_coef(list) -- Chebyshev polynomial coefficients that represent
            the stretch of the spectrum to put onto the fiducial index
            array (fix)
        lam_coeff(list) -- Chebyshev polynomial coefficients that convert
            fix to wavelength in nm

    Returns:
        Tuple containing the wavelength and spectrum
        (wavelength, spectrum) this is in the units of the fit coefficients.
        
    '''
    
    gix = scale_on_547(gridded)
    gsp = gridded.specw
    newix = chebval(fix, rgrd_coef)

    CC = chebfit(newix, fix, 3)
    XX = chebval(gix, CC)
    LL = chebval(XX, lam_coef)

    return (LL, gsp)



def plot_grid(grid, Xes=None):
    
    pl.figure(2) ; pl.clf()
    pl.figure(1) ; pl.clf()
    pl.ioff()
    ixs = []
    for index, g in enumerate(grid):
        if index < 800: continue
        if index > 900: continue
        if g.specw is None: continue
        if len(g.specw) == 0: continue
        if g.hgcoef is None: continue
        ix = np.arange(len(g.specw))
        ix_547 = np.argmin(np.abs(chebval(ix, g.hgcoef) - 547))

        if ix_547 > 180: continue
        if ix_547 < 100: continue

        ixs.append(ix_547)

        ix = np.arange(len(g.specw)) - ix_547

        add = 0
        if Xes is not None:
            try: add = Xes[index].specw
            except: pass
            if add is None: add = 0
        pl.figure(2)
        pl.plot(ix, g.specw+add, 'x')


        #if index > 500: break
    pl.show()

    pl.figure(3)
    pl.clf()
    pl.ion()
    pl.hist(ixs,50)


def rough_grid(extractions, lines=[365.0, 404.6, 435.8, 546.1, 578], 
    outname='rough_wavelength.npy'):
    '''Shift the extractions onto a coarse common pixel grid'''

    pl.figure(1)
    pl.clf()

    update_rate = len(extractions) / Bar.setup()
    #resfun = NPK.Fit.mpfit_residuals(NPK.Fit.sedm_wavelen)
    for ix,ex in enumerate(extractions):
        if not ex.ok: continue

        if ix % update_rate == 0: Bar.update()

        
        xs = []
        ys = []
        for line in lines:
            if line in ex.hg_lines:
                X = ex.hg_lines[line]
                if X + ex.xrange[0] > 2048: 
                    ex.ok = False
                    continue
                xs.append(ex.hg_lines[line])
                ys.append(line)

        xs = np.array(xs)
        ys = np.array(ys)

        if len(xs) == 0: 
            print ex.xrange, ex.yrange, ex.hg_lines
            continue
        
        coef = chebfit(xs,ys,2)
        vals = chebval(xs, coef)

        ex.hgcoef = coef
        err = (vals - ys)
        
        ix = np.arange(len(ex.specw))


    pl.ion()
    pl.xlim(360,600)
    pl.show()
    Bar.done()

    try:np.save(outname, extractions)
    except: pass
    return extractions



def RMS(vec):
    return np.sqrt(np.sum((vec-np.mean(vec))**2))

def fit_spectra_Hg_Xe(Hgs, Xes, kdtree, kdseg_ids, plot=False, outname='fit_spectra'):
    
    assert(len(Hgs) == len(Xes))

    if plot:
        pl.figure(1)
        pl.clf()
        pl.figure(2)
        pl.clf()
        pl.figure(1)

    ixs = np.arange(0,260,.01)
    hg_segids = [x.seg_id for x in Hgs]
    for spec_ix in xrange(len(Hgs)):
        hg = Hgs[spec_ix]
        xe = Xes[spec_ix]
        seg_id = hg.seg_id

        coeff = hg.hgcoef
        if coeff is None: 
            print spec_ix
            continue

        prev_rms = np.inf
        rms = np.inf
        for i in xrange(25):

            offsets = measure_offsets(hg.specw+xe.specw, coeff, plot=False)

            keys = np.array(offsets.keys())
            vals = np.array(offsets.values())
            #print np.sort(keys+vals)

            lams = chebval(ixs,coeff)
            pixs = []
            meas = []

            for lam, off in offsets.iteritems():
                if (lam+off) > 1000: continue
                if (lam+off) < 340: continue
                ix = ixs[np.argmin(np.abs(lam-lams))]
                pixs.append(ix)
                meas.append(lam+off)

            if np.any(np.isfinite(meas) == False): break
            if np.any(pixs == 0): break
            if len(pixs) < 2: break

            if i < 2: coeff = chebfit(pixs, meas, 3)
            else: coeff = chebfit(pixs, meas, 4)
            diff = chebval(pixs, coeff) - meas
            rms = RMS(diff)
            #print "%4.0i %6.3f" % (spec_ix, rms)

            if not np.isfinite(rms): pdb.set_trace()

            if (rms < 0.35): break
            if np.abs(rms - prev_rms) < 0.02: break
            prev_rms = rms


        if plot:
            pl.figure(2)
            lam = chebval(np.arange(len(hg.specw)), coeff)
            spec = hg.specw+xe.specw
            pl.plot(lam, spec)
            pl.figure(1)


        Hgs[spec_ix].lamcoeff = coeff
        Hgs[spec_ix].lamrms = rms
        print "-- %4.0i %6.3f" % (spec_ix, rms)
        #if rms > 3: pdb.set_trace()


    if plot:
        pl.show()
            

    np.save(outname, Hgs)
    return Hgs
        


def measure_offsets(spec, coeff, plot=False):
    '''Measure wavelength offsets of Hg and Xe extractions

    First uses a crude peak finding code to identify the locations of the 
    Xe lines.

    Then, the two Xe complexes are fit individually
        
    Returns:
        A dictionary of {line_wavelength_nm: offset_nm}.
    '''

    preffun = lambda x: ((x[1])**2 + (x[0]-120)**2)/1e4 + 1.
    resfun830 = NPK.Fit.mpfit_residuals(xe_830nm, preffun=preffun)

    # Preference function is designed to bias against adjusting the 
    preffun = lambda x: ( np.abs(x[0]-100)/100 + np.abs(x[1])/30 + np.abs(x[4]-100)/100)/3e4 + 1.0
    resfun890 = NPK.Fit.mpfit_residuals(xe_890nm, preffun=preffun)
    resfung = NPK.Fit.mpfit_residuals(NPK.Fit.gaussian5)

    #pk_pxs = SG.find_peaks_cwt(spec, np.arange(8,25))
    pk_pxs = SG.argrelmax(spec, order=10)[0]
    pk_lam = chebval(pk_pxs, coeff)
    ix = np.arange(len(spec))
    ll = chebval(ix, coeff)

    if plot:
        pl.figure(3)
        pl.clf()
        pl.plot(ll, spec)
        for p in pk_lam: 
            print "-", p
            pl.axvline(p,color='blue')


    offsets = {}

    # resfung parameters:

    for index, lam in enumerate([365.0, 404.6, 435.8, 546.1, 578]):
        roi = (ll > (lam-9)) & (ll < (lam+9))

        toll = ll[roi]
        tofit = spec[roi]
        if len(toll) == 0: continue
        offset = np.min(tofit)
        slope = 0
        sigma = 2.5
        scale = np.max(tofit)-offset

        p = [scale, lam, sigma, offset, slope]
        parinfo=[{'value': scale, 'limited':[1,0], 'limits':[0,0]}, 
            {'value': lam, 'limited':[1,1], 'limits':[lam-15,lam+15]}, 
            {'value': sigma, 'limited':[1,1], 'limits': [1.5,3.5]},
            {'value': offset, 'limited':[1,0], 'limits': [0,0]}, 
            {'value': slope}]
        res = NPK.Fit.mpfit_do(resfung, toll, tofit, parinfo)

        if res.status > 0:
            offsets[lam] = lam-res.params[1]


    # 830
    guess_pos = np.argmin(np.abs(pk_lam - 830))
    guess_lam = pk_lam[guess_pos]
    guess_lam = pk_lam[1]
    ok = ((guess_lam-35) < ll) & (ll < (guess_lam+45))

    if np.any(ok):
        sigma = 80
        lamoffset = guess_lam-830
        offset = np.min(spec[ok])
        peak = np.max(spec[ok]) - offset
        p = [sigma,lamoffset,offset,peak]

        parinfo = [{'value': sigma, 'limited':[1,1], 'limits':[25,250]}, 
                        {'value': lamoffset, 'limited':[1,1], 
                            'limits':[-50,50]},
                        {'value': offset,'limited':[1,0],'limits':[0,0]}, 
                        {'value': peak,'limited':[1,0],'limits':[0,0]}]

        res = NPK.Fit.mpfit_do(resfun830, ll[ok], spec[ok], parinfo)
        thell = ll[ok]
        ps = res.params.copy()
        ps[1] = 0
        thespec = xe_830nm(ps, thell)
        cl = np.sum(thespec*thell)/np.sum(thespec)
        offsets[cl] = -res.params[1]
        
        if plot:
            pl.plot(ll[ok], xe_830nm(res.params, ll[ok]),'x-')

    # Now 920
    guess_pos = np.argmin(np.abs(pk_lam - 920))
    guess_lam = pk_lam[guess_pos]
    guess_lam = pk_lam[0]
    ok = ((guess_lam-55) < ll) & (ll < (guess_lam+100))

    if np.any(ok):
        sigma = 100
        lamoffset = guess_lam-890
        if lamoffset > 5: lamoffset = 5
        if lamoffset < -5: lamoffset = -5
        offset = np.min(spec[ok])
        peak = np.max(spec[ok]) - offset
        p = [sigma,lamoffset,offset,peak, 500]
        p930 = 500

        parinfo = [     {'value': sigma, 'limited':[1,1], 'limits':[25,260]}, 
                        {'value': lamoffset, 'limited':[1,1], 'limits':[-50,50]},
                        {'value': offset,'limited':[1,0],'limits':[0,0]}, 
                        {'value': peak,'limited':[1,0],'limits':[0,0]},
                        {'value': p930,'limited':[1,0],'limits':[0,0]}]

        res = NPK.Fit.mpfit_do(resfun890, ll[ok], spec[ok], parinfo)

        thell = ll[ok]
        ps = res.params.copy()
        #ps[1] = 0
        thespec = xe_890nm(ps, thell)
        cl = np.sum(thespec*thell)/np.sum(thespec)
        offsets[cl] = -res.params[1]

        if plot:
            print res.status
            print res.niter, res.perror
            print res.params
            print res.debug
            print res.errmsg
            pl.plot(ll[ok], xe_890nm(res.params, ll[ok]),'x-')


    if plot: pdb.set_trace()
    return offsets




def read_spec_loc(spec_loc_fname):
    '''Returns the structure for the location of each spectrum'''


    return np.load(spec_loc_fname)

def xe_830nm(p, lam):
    '''Xe comoplex near 830 nm.

    See: http://www2.keck.hawaii.edu/inst/lris/arc_calibrations.html'''

    sigma,lamoffset,offset,peak = p

    sig=(2000*np.exp(-(lam-820.63-lamoffset)**2/sigma) + 
        1000*np.exp(-(lam-828.01-lamoffset)**2/sigma) + 
        100*np.exp(-(lam-834.68-lamoffset)**2/sigma) + 
        180*np.exp(-(lam-840.82-lamoffset)**2/sigma))
    sig = sig/max(sig)*peak + offset

    return sig


def save_fitted_ds9(fitted, outname='fine'):

    xs = []
    ys = []
    xvals  = np.arange(250)
    ds9  = 'physical\n'
    for ix,S in enumerate(fitted):
        if not S.__dict__.has_key('lamcoeff'): continue

        res = chebval(xvals, S.lamcoeff)
        invcoeffs = chebfit(res, xvals, 4)
        pxs= chebval([400., 500., 600., 656.3, 700., 800., 900.], invcoeffs)

        startx = S.xrange[0]


        for ix, px in enumerate(pxs):
            X = px + startx
            Y = np.poly1d(S.poly)(X)

            
            if ix == 3:
                ds9 += 'point(%s,%s) # point=cross text={%s}\n' % \
                    (X,Y, S.seg_id)
            else:
                ds9 += 'point(%s,%s) # point=cross\n' % (X,Y)


    f = open(outname+".reg", "w")
    f.write(ds9)
    f.close()
            
    





def xe_890nm(p, lam):
    '''Xe comoplex near 890 nm.

    See: http://www2.keck.hawaii.edu/inst/lris/arc_calibrations.html'''

    sigma,lamoffset,offset,peak,p937 = p

    sig=(5000*np.exp(-(lam-881.94-lamoffset)**2/sigma) + 
        1000*np.exp(-(lam-895.22-lamoffset)**2/sigma) + 
        1000*np.exp(-(lam-904.54-lamoffset)**2/sigma) + 
        1900*np.exp(-(lam-916.26-lamoffset)**2/sigma) +
        p937*np.exp(-(lam-937.42-lamoffset)**2/sigma) )
        #0000*np.exp(-(lam-951.33-lamoffset)**2/sigma) +
        #000*np.exp(-(lam-979.69-lamoffset)**2/sigma))
    sig = sig/max(sig)*peak + offset

    return sig


def assign_fit_to_spectra(target, gridded, rss, fix, stretch, lamfit):
    ''' Put wavelength solution into target
    
    Args:
        target(list of Extraciton): Target list of extractions
        gridded(list of Extractions[Ne]): Gridded spectra
        rss (list of float[Ne]): RSS wavelength fit error
        fix (list of float): Fiducial spectral index  for "stretch"
        stretch (list of float [Ne x Nfloat]): Float polynomials to stretch
        spectra onto
        lamfit (list of float [Ne x Ncoeff]): Wavelength fitting function
        coefficients
        
    Returns:
        New grid with coefficients assigned to it

    '''
    
    for i in xrange(len(gridded)):

        try: L,S = coeffs_to_spec(fix, gridded[i], stretch[i], lamfit[i])
        except: continue

        cc = chebfit(np.arange(*target[i].xrange), L, 5)
        target[i].lamcoeff = cc
        target[i].lamrms = rss[i]

    return target
   

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
    '''Wavelength.py performs:

        1. Rough wavelength calibration based on Mercury lamp lines.
        2. Fine wavelength calbiration based on Xenon and Mercury lamp lines.
        3. Extraction of spectra with a fine wavelength calibration.

        For each step, various report files are written, often ds9 region 
        files.

        the _coarse.npy file contains a length-1 list with a dictionary. The
            dictionary contains {wavelength: [(XY point)]} with all XY points
            for a given Hg emission line wavelength.

        the assoc_Hg.npy file contains a length-1 list with the associated
            Hg solutions from the above *_coarse file. 

        --dome comes from FindSpectra.py

        for rough set --dome [npy], --hgcat [txt], and --outname 
        for fine set --xefits [fits], --hefits [fits] --cdfits [fits] --hgassoc [npy], and --outname
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('step', type=str, help='One of [rough|fine|extract]')
    parser.add_argument('--infile', type=str, help='Infile name, purpose depends on step')
    parser.add_argument('--dome', type=str, help='Dome segment definition npy file. Used in rough.')
    parser.add_argument('--hgfits', type=str, help='Name of mercury fits file')
    parser.add_argument('--xefits', type=str, help='Name of xenon fits file')
    parser.add_argument('--cdfits', type=str, help='Name of cadmium fits file')
    parser.add_argument('--hefits', type=str, help='Name of helium fits file')
    parser.add_argument('--hgcat', type=str, help='Name of mercury sextractor catalog')
    parser.add_argument('--coarse', type=str, help='Name of coarse Hg solution [npy]')
    parser.add_argument('--hgassoc', type=str, help='Name of coarse Hg solution [npy]')
    parser.add_argument('--fine', type=str, help='Name of fine solution [npy]')
    parser.add_argument('--toextract', type=str, help='Name of fine solution [npy]')
    parser.add_argument('--outname', type=str, help='Prefix name of output file')

 
    args = parser.parse_args()
    outname = args.outname

    if args.step == 'rough':
        hgfits = args.hgfits
        catname = args.hgcat
        catalog = read_catalog(catname)
        spec_loc_fname = args.dome
        spec_loc = read_spec_loc(spec_loc_fname)
        hg_spec = find_hg_spectra(catalog, outname=outname)
        assoc_hg_with_flats(spec_loc, hg_spec)
    elif args.step == 'fine':
        XeDat = pf.open(args.xefits)
        HgDat = pf.open(args.hgfits)
        if args.cdfits is not None: CdDat = pf.open(args.cdfits)
        else: CdDat = None
        if args.hefits is not None: HeDat = pf.open(args.hefits)
        else: HeDat = None
        assoc_hg_spec = np.load(args.hgassoc)[0]

        Xe_E = extract(XeDat, assoc_hg_spec, filename="Xe_ext_"+args.outname)
        Hg_E = extract(HgDat, assoc_hg_spec, filename="Hg_ext_"+args.outname)

        Cd_E = None
        if CdDat is not None:
            Cd_E = extract(CdDat, assoc_hg_spec, filename="Cd_ext_"+args.outname)
        He_E = None
        if HeDat is not None:
            He_E = extract(HeDat, assoc_hg_spec, filename="He_ext_"+args.outname)

        gridded = rough_grid(Hg_E)

        stretchset, fiducial, xcors = stretch_set(gridded, Xe_E)
        fix, f1, f2 = fiducial
        fits, residuals, rss, locs = fit_all_lines(fix,
                                                    gridded,
                                                    Xe_E,
                                                    stretchset,
                                                    cd_spec = Cd_E,
                                                    he_spec = He_E)

        result = assign_fit_to_spectra(Hg_E, gridded, rss, 
                                        fix, stretchset, fits)

        result = median_fine_grid(result)
        print "Snapping solution into place"
        snap_solution_into_place_all(result, Hg_E, Xe_E, Cds=Cd_E, Hes=He_E)
        np.save(args.outname, result)


    elif args.step == 'fineold':
        ''' This step kept for historical purposes '''
        XeDat = pf.open(args.xefits)
        HgDat = pf.open(args.hgfits)

        assoc_hg_spec = np.load(args.hgassoc)[0]
        Xe_E = extract(XeDat, assoc_hg_spec, filename="Xe_ext_"+args.outname)
        Hg_E = extract(HgDat, assoc_hg_spec, filename="Hg_ext_"+args.outname)
        #Xe_E = np.load('Xe_ext_tester.npy')
        #Hg_E = np.load('Hg_ext_tester.npy')
        gridded = rough_grid(Hg_E, outname=outname)
        import pdb
        pdb.set_trace()

        tree, segments = hg_to_kdtree(assoc_hg_spec)
        fitted = fit_spectra_Hg_Xe(gridded, Xe_E, tree, segments,
            outname=outname)
        fitted = median_fine_grid(fitted)
        np.save("%s.npy" % outname, fitted)
        save_fitted_ds9(fitted, outname=outname)
    elif args.step == 'extract':
        
        fitted = np.load(args.fine)
        hdu = pf.open(args.toextract)
        outname = args.outname

        ww = wavelength_extract(hdu, fitted, filename=outname)
        
    sys.exit()
    
    
    #hg_spec = np.load('Hg.txt.npy')[0]
    ##hg_spec = find_hg_spectra(catalog, outname=outname)
    ##assoc_hg_spec = assoc_hg_with_flats(spec_loc, hg_spec)
    #assoc_hg_spec = np.load('assoc_Hg.npy')[0]

    #XeDat = pf.open("Xe.fits")
    ##S2 = extract(XeDat, assoc_hg_spec, filename="raw_xe_extractions.npy")
    #S2 = np.load("raw_xe_extractions.npy")
    ##S = extract(dat, assoc_hg_spec)
    #extractions = np.load('raw_hg_extractions.npy')

    ##gridded = rough_grid(extractions)
    #gridded = np.load('gridded.npy')
    ##fitted = fit_spectra_Hg_Xe(gridded, S2, plot=False)
    #fitted = np.load('fitted.npy')

    #save_fitted_ds9(fitted)


