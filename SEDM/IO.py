import numpy as np
import scipy as sp
import scipy.io
import astropy
from astropy.table import Table as T


def load_cube(fname):
    
    cub = sp.io.loadmat(fname)

    return cub['Cube'][0][0]

def load_spec(fname):

    spec = sp.io.loadmat(fname)
    return spec["Spec"][0][0]

def write_spec(spec, fname):

    t = T([spec["Wave"][:,0], spec["Source"][:,0], spec["BS_Source"][:,0]],
        names=("Wave", "Source", "BS Source"))

    t.write(fname, format="ascii.commented_header")
