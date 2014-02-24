# Code by Nick Konidaris
# (c) 2014
#
# nick.konidaris@gmail.com
# For the SED Machine Spectrograph
#

import IO
import numpy as np
import scipy as sp
import pylab as pl

import scipy.signal
import scipy.spatial
from scipy.interpolate import interp1d

from SEDM import Flexure as FF

reload (IO)

def lambda_to_dlambda(lam):
    """Extraploate dlambda of array"""

    dlam = np.diff(lam)
    ok = np.isfinite(dlam)
    xx = np.arange(len(dlam))
    try:
        dlamf = np.poly1d(np.polyfit(xx[ok],
                                dlam[ok],
                                4))
    except:
        import pdb
        pdb.set_trace()

    dlam = dlamf( np.arange(len(lam)))

    return dlam



def segmap_to_kdtree(SegMap):
    """Creates a segmentation map KD Tree

    Args:
        SegMap: The segmentation map

    Returns a scipy.spatial.KDTree
    """
    data = []
    oks = []
    for i in xrange(len(SegMap)):
        seg = SegMap[i][0]

        x = seg['OnSkyX'][0][0]
        y = seg['OnSkyY'][0][0]

        if np.isfinite(x) and np.isfinite(y):
            oks.append(i)
            data.append( (x,y ))


    kt = scipy.spatial.KDTree(np.array(data))
    return kt, oks

def spectra_in_annulus(KT, SegMap, X, Y, small=200, large=300, ixmap=None,
                        pixel_shift=0):
    """Queries the KD tree for points within the annulus, returns the spectra.

    Args:
        KT: The KD tree. This can be created in segmap_to_kdtree
        SegMap: The segmentation map
        X,Y: The X/Y position of the spectrum to extract
        small,large: The near and far distances of the annulus
        ixmap: The map of KT.data to SegMap
        pixel_shift: For flexure correction, the number of pixels to shift
            the wavelength solution

    Returns a list of (Wavelength, Spectra). E.g.:
        [   [array(365 .. 1000) , array( 50 .. 63)] ,
            [array(366 .. 1001) , array( 53 .. 61)] ,
            ]

    """
    near = set(KT.query_ball_point( (X, Y), small))
    far = set(KT.query_ball_point( (X, Y), large))

    if ixmap is None: ixmap = range(len(SegMap))

    in_annulus = far.difference(near)
    results = []
    for index in in_annulus:
        ix = ixmap[index]
        lam = SegMap[ix]['WaveCalib'][0][0][:]
        lam = FF.shift_pixels(lam, pixel_shift)
        results.append( (lam,
                        SegMap[ix]['SpexSpecFit'][0][0],
                        ix))

    return results



def sky_median(KT, SegMap, X=2, Y=10, distance=7, ixmap=None, pixel_shift=0):
    """Generates a median spectrum estimating the "sky".

    Args:
        KT: The KD tree. This can be created in segmap_to_kdtree
        SegMap: The segmentation map
        X,Y: The X/Y position of the spectrum to extract. Default is middle.
        distance: Distance to the spaxel in pixels. Default is 900
        pixel_shift: For flexure correction, the number of pixels to shift
            the wavelength solution

    Returns:
        {'wave_nm': wavelength of sky spectrum, 
            'spec_adu', the median sky spectrum, 
            'all_spec': and a matrix of spectra,
            'num_spec': The number of spectra}

    """

    specs = spectra_near_position(KT, SegMap, X, Y, distance, ixmap=None,
        pixel_shift=pixel_shift)
    return interp_and_sum_spectra(specs)

def spectra_near_position(KT, SegMap, X, Y, distance=2, ixmap=None,
                            pixel_shift = 0):
    """Queries the KD tree for points near the source, returns the spectra.

    Args:
        KT: The KD tree. This can be created in segmap_to_kdtree
        SegMap: The segmentation map
        X,Y: The X/Y position of the spectrum to extract
        distance: Distance to the spaxel in pixels. Default is 100
        ixmap: The mapping of index from KT.data to SegMap
        pixel_shift: For flexure correction, the number of pixels to shift
            the wavelength solution

    Returns a list of (Wavelength, Spectra, SegMap index). E.g.:
        [   [array(365 .. 1000) , array( 50 .. 63)] ,
            [array(366 .. 1001) , array( 53 .. 61)] ,
            ]

    """

    ixs = KT.query_ball_point( (X, Y), distance)

    if ixmap is None:
        ixmap = range(len(SegMap))

    results = []
    for index in ixs:
        ix = ixmap[index]
        if len(SegMap[ix]['WaveCalib'][0]) == 0:
            continue
        lam = SegMap[ix]['WaveCalib'][0][0][:]
        lam = FF.shift_pixels(lam, pixel_shift)

        results.append( (lam,
                        SegMap[ix]['SpexSpecFit'][0][0], 
                        ix))

    return results

def interp_and_sum_spectra(specs, sky_spec = None, onto=None):
    """Sums and interpolates spectra to a common grid.

    The result of spectra_near_pixel are interpolated to a common grid and summed.
    If a sky spectrum is provided it is subtracted off of each spectrum.

    Args:
        specs: List of spectra. Each spectrum is a (wavelength, flux) tuple.
        sky_spec: (Wave, Spec) of sky spectrum to be subtracted off
        onto: The wavelengths on which to interpolate the spectrum

    Returns:
         {'wave_nm': wavelength, # in nm
            'spec_adu': spectrum, # in adu
            'num_spec': total # of spectra participating
            'all_spec': array of [nwave, nspec] interpolated onto same grid
            'segments': array of [nspec] segment numbers}
    """

    # for now, assume the first spectrum has a good wavelength range
    # TODO: challenge assumption in comment above

    wave = specs[0][0][::-1]
    spec = specs[0][1][::-1]
    minw = np.nanmin(wave)
    maxw = np.nanmax(wave)

    if onto is not None:
        # Interpolate onto a specified grid
        specfun = interp1d(wave, spec, bounds_error=False, fill_value=0)
        spec = specfun(onto)
        wave = onto[:]

    all = np.zeros((len(wave) , len(specs)))

    if sky_spec is not None:
        sky_wave, sky_spec = sky_spec
        dlamsky = lambda_to_dlambda(sky_wave)
        skyf = interp1d(sky_wave, sky_spec/dlamsky, 
            bounds_error=False, fill_value=0.0)

    nspec = np.ones(len(spec), dtype=np.int)
    for i in xrange(1, len(specs)):
        lam = specs[i][0][::-1]
        ss = specs[i][1][::-1]

        dlam = lambda_to_dlambda(lam)
        dlamf = interp1d(lam, dlam, bounds_error=False,
            fill_value=0.0)

        ok = np.isfinite(lam) & np.isfinite(ss) & \
            (lam > minw) & (lam < maxw)
        if not np.any(ok):
            print "passing"
            continue

        interp_fun = interp1d(lam[ok], ss[ok], bounds_error=False,
            fill_value = np.nan)

        if sky_spec is None:
            the_spec = interp_fun(wave)
        else:
            the_spec = interp_fun(wave) - skyf(wave)*dlamf(wave)

        all[:,i] = the_spec

        ok2 = np.isfinite(the_spec)
        nspec[ok2] += 1

    return  {'wave_nm': wave,
        'spec_adu': np.nansum(all, axis=1),
        'num_spec': nspec,
        'all_spec': all}


def segmap_to_img(SegMap, sky_spec=None, minl=400, maxl=850):
    """Take seg.Mean[XY] and return an image of the Segmentation Map.
    
    Creates a tesselated hexagon image of the segmentation map. 

    Args:
        seg: The segmentation map
        sky_spec: A 2-tuple of (Wavelength, Spectrum) representing the night
            sky spectrum per spaxel. If none does nothing, if set subtracts
            from the spectrum before summing.   
        minl/maxl: The minimum and maximum wavelength over which to take
            the median

    Returns:
        A three-tuple of numpy array. Format is (X positions, Y positions, Values)
        Values represent median value of the spectrum between minl and maxl
    
    """

    Xs = []
    Ys = []
    Values = []
    if sky_spec is not None:
        sky_wave, sky_spec = sky_spec
        dlamsky = lambda_to_dlambda(sky_wave)

        skyf = interp1d(sky_wave, sky_spec/dlamsky, bounds_error=False,
            fill_value=0.0)

        dlamskyf = interp1d(sky_wave, dlamsky, bounds_error=False,
            fill_value=0.0)

        minw = np.nanmin(sky_wave)
        maxw = np.nanmax(sky_wave)
        
    for i in xrange(len(SegMap)):
        seg = SegMap[i][0]
    
        if len(seg['WaveCalib']) == 0: continue

        wave = seg['WaveCalib'][0][::-1]
        spec = seg['SpexSpecFit'][0][::-1]

        if len(wave) < 10: continue

        roi = (wave > minl) & (wave < maxl)

        Xs.append(seg['OnSkyX'][0][0])
        Ys.append(seg['OnSkyY'][0][0])

        if sky_spec is None:
            Values.append(np.median(spec[roi]))
        else:
            try:
                ok = np.isfinite(wave) & np.isfinite(spec) & \
                    (wave > minw) & (wave < maxw) & (wave > minl) & \
                    (wave < maxl)
            except:
                import pdb
                pdb.set_trace()
            if not np.any(ok):
                Values.append(np.median(skyf(wave)))
                continue
            
            dlam = lambda_to_dlambda(wave)
            local_sky = skyf(wave) * dlamskyf(wave)

            #local_sky = skyf(wave[ok])
            Values.append(np.median(spec[ok] - local_sky[ok]))

    return map(np.array, (Xs, Ys, Values))





def cub_to_img(cub):
    shp = cub["Cube"].shape
    skys = med_spec(cub)

    img = np.zeros((shp[1], shp[2]))

    for i in xrange(shp[1]):
        for j in xrange(shp[2]):
            ms = sp.signal.medfilt(cub["Cube"][:,i,j], 3)
            ms = cub["Cube"][:,i,j]
            img[i,j] = np.sum(ms - skys)

    return img

def go(fname):
    cub = IO.load_cube(fname)

    img = cub_to_img(cub)

    return img 

def med_spec(cub):
    shp = cub["Cube"].shape
    AllSpec = np.zeros((shp[0], shp[1] * shp[2]))

    cnt = 0
    for i in xrange(shp[1]):
        for j in xrange(shp[2]):
            s = cub["Cube"][:,i,j]
            if np.any(np.isfinite(s)):
                AllSpec[:, cnt] = cub["Cube"][:,i,j]
                cnt += 1
            
    ms = np.median(AllSpec[:,0:cnt],axis=1)
    print cnt
    return ms
    

def avg_spec(cub):

    shp = cub["Cube"].shape
    Spec = np.zeros(shp[0])
    Nok = np.zeros(shp[0])

    for i in xrange(shp[1]):
        for j in xrange(shp[2]):
            s = cub["Cube"][:,i,j]
            ok = np.where(np.isfinite(s))
            Spec[ok] += s[ok]
            Nok[ok] += 1

    pl.plot(Spec/Nok)

    return Spec/Nok


def extract(cub, x, y, aper=[5,8,12]):

    x,y = map(np.round, [x,y])

    shp = cub["Cube"].shape

    XX, YY = np.meshgrid(np.arange(shp[1]), np.arange(shp[2]))


    d = np.sqrt(((XX-x)**2 + (YY-y)**2))

    obj = np.where(d <= aper[0])
    sky = np.where((aper[1] <= d) & (d <= aper[2]))

    
    Object = np.zeros(shp[0])
    Sky = np.zeros(shp[0])
    Nobj = len(obj[0])
    Nsky = len(sky[0])

    for i in xrange(Nobj):
        x,y = obj[0][i],obj[1][i]
        Object += cub["Cube"][:,x,y]

    for i in xrange(Nsky):
        x,y = sky[0][i],sky[1][i]
        Sky += cub["Cube"][:,x,y]

    print Nsky, Nobj

    pl.figure(3)
    pl.clf()
    pl.plot(Object - Sky/Nsky*Nobj)
    

