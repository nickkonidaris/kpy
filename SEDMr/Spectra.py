


import datetime
import numpy as np
import json
import pyfits
import scipy.io
import matplotlib.pyplot as pl
from matplotlib.backend_bases import KeyEvent 
from matplotlib.backend_bases import PickEvent
import scipy, scipy.spatial
from numpy.polynomial.chebyshev import chebfit, chebval

def find_ha(cc):
    ix = np.arange(30, 100, .1)

    haix = np.argmin(np.abs(chebval(ix, cc) - 656.3))

    return ix[haix]

class Spectra(object):
    KT = None # The KD Tree
    data = None # The actual data
    good_positions = [] # The mapping of KT data to data so that
                        # some ix in the KT.data is the same
                        # as data[KT.good_positions[ix]]

    def __init__(self, data=None, minl=500, maxl=700):
        
        positions = []
        good_positions = []
        self.data = data

        for ix, el in enumerate(data):
            
            try:
                l,fl = el.get_counts()
                fl = el.specw
                haix = find_ha(el.lamcoeff)
            except: continue
            
            good_positions.append(ix)

            x = el.X_as
            y = el.Y_as

            positions.append( (x,y) )

        if len(positions) == 0:
            raise Exception("For some reason, no good spectrum exists in the submitted spectral set.")

        data = np.array(positions)
        bad = (data != data)
        data[bad] = -999
        self.KT = scipy.spatial.KDTree(data)
        self.good_positions = np.array(good_positions)

    def to_xyv(self, lmin=500, lmax=700, coefficients='lamcoeff'):
        ''' Method convers a spetral set into X, Y, and value positions
        
        X is the location of Ha based on either lambda coefficients or
            median coefficients
        Y is based on the segmentation map.
        The value is the median signal strength in the lmin to lmax range.

        Returns:
            X,Y,V tuple representing the X location (in arcsec), Y location
                (in arcsec), and the median value (V) of the spaxel
                between lmin and lmax.
        '''

        allowed_coeff = ['lamcoeff', 'mdn_coeff', 'hgcoef']
        if coefficients not in allowed_coeff:
            raise Exception("Named coefficient %s should be one of %s" % 
                (coefficients, allowed_coeff))

        
        Xs = []
        Ys = []
        Vs = []

        for ix, XY in enumerate(self.KT.data):

            datix = self.good_positions[ix]
            el = self.data[datix]
            try:
                l,fl = el.get_flambda()
                lamcoeff = getattr(el, coefficients)
                haix = find_ha(lamcoeff)
            except: 
                continue

            ok = (l > lmin) & (l <= lmax)

            Xs.append(XY[0])
            Ys.append(XY[1])
            Vs.append(np.median(el.spec[ok]))
        

        return (np.array(Xs),
                np.array(Ys),
                np.array(Vs))
