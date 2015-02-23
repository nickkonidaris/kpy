


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
                l,fl = el.get_flambda()
                haix = find_ha(el.lamcoeff)
            except: continue
            
            good_positions.append(ix)

            x = el.X_as
            y = el.Y_as

            positions.append( (x,y) )

        data = np.array(positions)
        self.KT = scipy.spatial.KDTree(data)
        self.good_positions = np.array(good_positions)

    def to_xyv(self, lmin=500, lmax=700):
        
        Xs = []
        Ys = []
        Vs = []

        for ix, XY in enumerate(self.KT.data):

            datix = self.good_positions[ix]
            el = self.data[datix]
            try:
                l,fl = el.get_flambda()
                haix = find_ha(el.lamcoeff)
            except: continue

            ok = (l > lmin) & (l <= lmax)

            Xs.append(XY[0])
            Ys.append(XY[1])
            Vs.append(np.median(el.spec[ok]))
        

        return (np.array(Xs),
                np.array(Ys),
                np.array(Vs))
