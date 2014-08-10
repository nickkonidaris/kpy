
import pdb
import numpy as np
import pylab as pl
import pyfits as pf
import sys


import NPK.Fit 
import NPK.Bar as Bar
from astropy.table import Table 

from scipy.spatial import KDTree 
import scipy.signal as SG


from numpy.polynomial.chebyshev import chebfit, chebval

import Extraction
reload(NPK.Fit)
reload(Extraction)

def read_catalog(catname):
    '''Read a sextractor catalog called catname and return'''


    cat = Table.read(catname, format='ascii.sextractor')

    return cat



def assoc_hg_with_flats(domedat, hgcat, guess_offset= {365.0: 231,
    404.65:214, 435.83:193, 546.07: 133, 578.00: 110}, outname='assoc_Hg'):
    
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
        tracefun = np.poly1d(spec_pos['poly'])
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

        
def find_hg_spectra(lines, dYlimit=3, outname="find_spectra"):
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
        results = kdt.query_ball_point(point, 33)
        X, Y = point
        num = line['NUMBER']
        flux = float(line['FLUX_ISO'])

        for resix in results:
            resX, resY = data[resix]
            resFlux = float(lines[resix]['FLUX_ISO'])
            dist = np.abs(resX-X)
            dX = resX - X
            
            dY = np.abs(resY - Y)
            if dY > dYlimit: continue

            ratio = flux/resFlux
            if (2 < ratio < 7) and (-16 < dX < -12):
                reg += 'circle(%s,%s, 3) # color=red\n' % (X, Y)
                reg += 'circle(%s,%s, 3) # color=magenta\n' % (resX, resY)
                hg546.append((X,Y))
                hg578.append((resX,resY))

            if (1/6. < ratio < 1) and (-24 < dX < -18):
                hg405.append((X,Y))
                hg436.append((resX,resY))
                reg += 'circle(%s,%s, 3) # color=green\n' % (X, Y)
                reg += 'circle(%s,%s, 3) # color=yellow\n' % (resX, resY)

            if (1 < ratio < 15) and (22.5 < dX < 33):
                hg365.append((resX,resY))
                reg += 'circle(%s,%s, 3) # color=blue\n' % (resX, resY)

    f = open(outname + ".reg", "w")
    f.write(reg)
    f.close()
    
    res = [{365.0: np.array(hg365),
            404.65: np.array(hg405),
            435.83: np.array(hg436),
            546.07: np.array(hg546),
            578.00: np.array(hg578)}]


    np.save(outname, res)

    return res[0]

def extract(HDUlist, assoc_hg_spec, filename='raw_hg_extractions'):


    dat = HDUlist[0].data

    
    extractions = []
    update_rate = len(assoc_hg_spec) / Bar.setup()

    for ix, ss in enumerate(assoc_hg_spec):
        if ix % update_rate == 0: Bar.update()
        if not ss['ok']: 
            extractions.append(Extraction.Extraction(seg_id=ss['seg_cnt'], 
                ok=False))
            continue
        minx = np.max((0,ss['xs'][0]-5))
        maxx = np.min((minx + 265,2047))
        yfun = np.poly1d(ss['poly'])

        xpos = xrange(minx, maxx)
        res = np.zeros(len(xpos))
        res[:] = np.nan
        for i in xrange(len(xpos)):
            X = xpos[i]
            Y = yfun(X)
            if not np.isfinite(X) or not np.isfinite(Y):
                continue
            Ys = slice(np.max((0,np.int(Y)-2)),np.min((np.int(Y)+2, 2047)))
            res[i] = np.sum(dat[Ys,X])

        hg_lines = {}
        for wave, pix in ss.iteritems():
            if type(wave) is float:
                hg_lines[wave] = pix-minx

        extractions.append( 
            Extraction.Extraction(xrange=(minx,maxx), yrange=(yfun(xpos[0]), 
                                    yfun(xpos[-1])),
                                    poly=ss['poly'], spec=res, 
                                    hg_lines=hg_lines, seg_id=ss['seg_cnt'],
                                    ok=True))

    Bar.done()

    np.save(filename, extractions)
    return extractions

def rough_grid(extractions, lines=[365.0, 404.65, 435.83, 546.07, 578.0], 
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
                if X + ex.xrange[0] > 2048: continue
                xs.append(ex.hg_lines[line])
                ys.append(line)

        xs = np.array(xs)
        ys = np.array(ys)

        if len(xs) == 0: continue
        
        coef = chebfit(xs,ys,3)
        vals = chebval(xs, coef)

        ex.hgcoef = coef
        err = vals - ys
        
        ix = np.arange(len(ex.spec))
        if len(xs) > 4: 
            pl.plot(chebval(ix, coef),ex.spec)


    pl.ion()
    pl.show()
    Bar.done()

    np.save(outname, extractions)
    return extractions

def RMS(vec):
    return np.sqrt(np.sum((vec-np.mean(vec))**2))

def fit_spectra_Hg_Xe(Hgs, Xes, plot=False):
    
    assert(len(Hgs) == len(Xes))

    if plot:
        pl.figure(1)
        pl.clf()
        pl.figure(2)
        pl.clf()
        pl.figure(1)

    ixs = np.arange(0,260,.01)
    for spec_ix in xrange(len(Hgs)):
        hg = Hgs[spec_ix]
        xe = Xes[spec_ix]
        if 'hgcoef' not in hg.__dict__: continue
        coeff = hg.hgcoef
        
        prev_rms = np.inf
        rms = np.inf
        for i in xrange(25):

            offsets = measure_offsets(hg.spec+xe.spec, coeff, plot=False)
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
            print "%4.0i %6.3f" % (spec_ix, rms)

            if not np.isfinite(rms): pdb.set_trace()

            if (rms < 0.35): break
            if np.abs(rms - prev_rms) < 0.02: break
            prev_rms = rms


        if plot:
            pl.figure(2)
            lam = chebval(np.arange(len(hg.spec)), coeff)
            spec = hg.spec+xe.spec
            pl.plot(lam, spec)
            pl.figure(1)

        Hgs[spec_ix].lamcoeff = coeff
        Hgs[spec_ix].lamrms = rms
        print "-- %4.0i %6.3f" % (spec_ix, rms)
        #if rms > 3: pdb.set_trace()


    if plot:
        pl.show()
            

    return Hgs
        


def measure_offsets(spec, coeff, plot=False):
    '''Measure wavelength offsets of Hg and Xe extractions

    First uses a crude peak finding code to identify the locations of the 
    Xe lines.

    Then, the two Xe complexes are fit individually
        
    Returns:
        A dictionary of {line_wavelength_nm: offset_nm}.
    '''
    resfun830 = NPK.Fit.mpfit_residuals(xe_830nm)
    resfun890 = NPK.Fit.mpfit_residuals(xe_890nm)
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

    for index, lam in enumerate([365.0, 404.65, 435.83, 546.07, 577.95]):
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

        parinfo = [{'value': sigma, 'limited':[1,1], 'limits':[25,150]}, 
                        {'value': lamoffset, 'limited':[1,1], 
                            'limits':[-30,30]},
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
        sigma = 80
        lamoffset = guess_lam-890
        offset = np.min(spec[ok])
        peak = np.max(spec[ok]) - offset
        p = [sigma,lamoffset,offset,peak]

        parinfo = [{'value': sigma, 'limited':[1,1], 'limits':[25,150]}, 
                        {'value': lamoffset, 'limited':[1,1], 
                            'limits':[-30,30]},
                        {'value': offset,'limited':[1,0],'limits':[0,0]}, 
                        {'value': peak,'limited':[1,0],'limits':[0,0]}]

        res = NPK.Fit.mpfit_do(resfun890, ll[ok], spec[ok], parinfo)
        thell = ll[ok]
        ps = res.params.copy()
        ps[1] = 0
        thespec = xe_890nm(ps, thell)
        cl = np.sum(thespec*thell)/np.sum(thespec)
        offsets[cl] = -res.params[1]

        if plot:
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


def save_fitted_ds9(fitted, outname='fitted.reg'):

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


    f = open(outname, "w")
    f.write(ds9)
    f.close()
            
    





def xe_890nm(p, lam):
    '''Xe comoplex near 890 nm.

    See: http://www2.keck.hawaii.edu/inst/lris/arc_calibrations.html'''

    sigma,lamoffset,offset,peak = p

    sig=(5000*np.exp(-(lam-881.94-lamoffset)**2/sigma) + 
        1000*np.exp(-(lam-895.22-lamoffset)**2/sigma) + 
        1000*np.exp(-(lam-904.54-lamoffset)**2/sigma) + 
        1900*np.exp(-(lam-916.26-lamoffset)**2/sigma) +
        1600*np.exp(-(lam-937.42-lamoffset)**2/sigma) )
        #0000*np.exp(-(lam-951.33-lamoffset)**2/sigma) +
        #000*np.exp(-(lam-979.69-lamoffset)**2/sigma))
    sig = sig/max(sig)*peak + offset

    return sig

    

if __name__ == '__main__':
    
    
    dat_fname = sys.argv[1]
    catname = sys.argv[2]
    spec_loc_fname = sys.argv[3]
    outname = sys.argv[4]

    dat = pf.open(dat_fname)

    #catalog = read_catalog(catname)
    spec_loc = read_spec_loc(spec_loc_fname)
    hg_spec = np.load('Hg.txt.npy')[0]
    #hg_spec = find_hg_spectra(catalog, outname=outname)
    #assoc_hg_spec = assoc_hg_with_flats(spec_loc, hg_spec)
    assoc_hg_spec = np.load('assoc_Hg.npy')[0]

    XeDat = pf.open("Xe.fits")
    #S2 = extract(XeDat, assoc_hg_spec, filename="raw_xe_extractions.npy")
    S2 = np.load("raw_xe_extractions.npy")
    #S = extract(dat, assoc_hg_spec)
    extractions = np.load('raw_hg_extractions.npy')

    #gridded = rough_grid(extractions)
    gridded = np.load('gridded.npy')
    #fitted = fit_spectra_Hg_Xe(gridded, S2, plot=False)
    fitted = np.load('fitted.npy')

    save_fitted_ds9(fitted)


