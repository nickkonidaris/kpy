
import argparse, os, pdb, sys
import numpy as np
import pylab as pl
import pyfits as pf
import scipy.signal as SG
from scipy.spatial import KDTree 

from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d
import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
import SEDMr.Spectra as SS

import sys
sys.setrecursionlimit(10000)


scale = 1.0
H2P = np.array([[np.sqrt(3), np.sqrt(3)/2], [0, 3/2.]]) * scale
P2H = np.array([[np.sqrt(3)/3, -1/3.], [0, 2/3.]]) / scale

# Rotation matrix
theta = np.deg2rad(-20)
ROT = np.array([[np.cos(theta), -np.sin(theta)], 
                [np.sin(theta),  np.cos(theta)]])


'''

Axial coordinates of a cube



      (q,r-1)   (q+1,r-1)

 (q-1, r)   (q, r)   (q+1, r)

      (q-1,r+1)   (q, r+1)
                

Even Q coordinates look like

    
        0,0     1,0     2,0     3,0
    0,1     1,1     2,1     3,1     4,1
        0,2     1,2     2,2     3,2


'''

def extraction_to_cube(exts, outname="G.npy"):
    ''' Convert the extraction to sky coordinates


    Args:
        exts: The list of extractions

    Results:
        outname: The file created with the new extraciton

    Returns:
        The new data cube with the following coordinate positions populated:

        X_as: X position in arcsecond
        Y_as: Y position in arcsecond
        Z_as: Z position in arcsecond (the Z coordinate is runs 45 degree to X
            and is not a ~3rd~ dimension).

        Q_as: The axial Q coordinate
        R_as: The axial R coordinate

        The relationship of Q/R to X/Y is defined through the pixel mapping
        matrix times the rotation matrix
            P2H = np.array([[np.sqrt(3)/3, -1/3.], [0, 2/3.]]) 
            Rot (22 degree)
    '''

    
    Xs = [None] * len(exts)
    Ys = [None] * len(exts)
    segids = [el.seg_id for el in exts]

    for ext in exts:
        ext.Q_as = None
        ext.R_as = None
    
    for ix, ext in enumerate(exts):
        Xs[ix] = -999
        Ys[ix] = -999

        if ext.lamcoeff is None: continue
        ixs = np.arange(*ext.xrange)
        LL = chebval(ixs, ext.lamcoeff)
        
        ix_ha = np.argmin(np.abs(LL-656.3))

        Xs[ix] = ixs[ix_ha]
        Ys[ix] = np.mean(ext.yrange)

    dat = np.array([Xs,Ys],dtype=np.float).T

    tree = KDTree(dat)

    ignore, Center = tree.query([1024,1024], 1)
    exts[Center].Q_as = 0
    exts[Center].R_as = 0


    def populate_hex(to_populate):
        ''' Breadth-first search 

            For each spaxel in the datacube this piece of code identified
            the relative offset based on the rotation matrix defined
            earlier in the file
        '''
        # Query 7 for the object + six surrounding members
        v = np.array([Xs[to_populate], Ys[to_populate]])
        Dists, Ixs = tree.query(v, 7)
        Tfm = P2H * ROT / np.median(Dists) * 1.8

        if Dists[0] < 2: 
            Dists = Dists[1:]
            Ixs = Ixs[1:]

        ok = Dists < 80
        Dists = Dists[ok]
        Ixs = Ixs[ok]

        q_this = exts[to_populate].Q_as
        r_this = exts[to_populate].R_as

        pl.figure(1)
        pl.clf()
        for nix in Ixs:
            nv = np.array([Xs[nix], Ys[nix]])
            D = nv-v
            rnd = np.round(np.dot(Tfm , D))

            if exts[nix].Q_as is None:
                exts[nix].Q_as = q_this + rnd[0]
                exts[nix].R_as = r_this + rnd[1]
                populate_hex(nix)

            pl.plot(D[0], D[1], 'o')
            pl.text(D[0], D[1], "%s" % rnd)

    populate_hex(Center)

    # Now convert Q/R to even-Q X/Y
    #

    Qs = np.array([ext.Q_as for ext in exts], dtype=np.float)
    Rs = np.array([ext.R_as for ext in exts], dtype=np.float)

    Xs = Qs
    Zs = Rs - (Qs + (Qs%2))/2.0
    Ys = -Xs - Zs

    for ix, ext in enumerate(exts):
        ext.X_as = Xs[ix]
        ext.Y_as = Ys[ix]
        ext.Z_as = Zs[ix]
        

    np.save(outname, exts)



parser = argparse.ArgumentParser(description=\
    '''Cube.py:

    Convert an extracted file into a data cube
        
    ''', formatter_class=argparse.RawTextHelpFormatter)


parser.add_argument('extracted', type=str, help='Extracted file')
parser.add_argument('--outname', type=str, help='Output cube name')

args = parser.parse_args()

if __name__ == '__main__':
    
    if args.outname is not None:
        args.outname = args.outname.rstrip('.npy')

    infile = args.extracted

    ext = np.load(infile)

    cube = extraction_to_cube(ext)





