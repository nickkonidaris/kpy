
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
theta = np.deg2rad(-37+13.5)
ROT = np.array([[np.cos(theta), -np.sin(theta)], 
                [np.sin(theta),  np.cos(theta)]])


'''

See figures here:


https://www.evernote.com/l/AA1s_YNLWMJBQrCCMZUn8-POndzDsbf4SwA

Axial coordinates of a cube



      (q,r-1)   (q+1,r-1)

 (q-1, r)   (q, r)   (q+1, r)

      (q-1,r+1)   (q, r+1)
                

Even Q coordinates look like

    
        0,0     1,0     2,0     3,0
    0,1     1,1     2,1     3,1     4,1
        0,2     1,2     2,2     3,2


'''

def QR_to_img(exts, Size=2, outname="cube.fits"):
    
    Qs = np.array([ext.Q_ix for ext in exts], dtype=np.float)
    Rs = np.array([ext.R_ix for ext in exts], dtype=np.float)

    minx = Size * np.sqrt(3) * (np.nanmin(Qs) + np.nanmin(Rs)/2)
    miny = Size * 3./2. * (np.nanmin(Rs))
    maxx = Size * np.sqrt(3) * (np.nanmax(Qs) + np.nanmax(Rs)/2)
    maxy = Size * 3./2. * (np.nanmax(Rs))


    Dx = maxx-minx
    Dy = maxy-miny
    l_grid = Wavelength.fiducial_spectrum()
    l_grid = l_grid[::-1]
    dl_grid = np.diff(l_grid)
    l_grid = l_grid[1:]

    img = np.zeros((Dx, Dy, len(l_grid)/2))
    img[:] = np.nan

    XSz = img.shape[0]/2
    YSz = img.shape[1]/2

    allspec = np.zeros((len(exts), len(l_grid)/2))
    for cnt, ext in enumerate(exts):
        if ext.xrange is None: continue
        if ext.exptime is None: ext.exptime = 1
        if ext.lamcoeff is None: continue

        ix = np.arange(*ext.xrange)
        l = chebval(ix, ext.lamcoeff)
        s = ext.specw

        f = interp1d(l, s, fill_value=np.nan, bounds_error=False)
        fi = f(l_grid) / dl_grid
        fi = fi[0:len(fi):2] + fi[1:len(fi):2]

        allspec[cnt, :] = fi

        try:
            x = np.round(Size*np.sqrt(3.) * (ext.Q_ix + ext.R_ix/2)) + XSz
            y = np.round(Size*3./2. * ext.R_ix) + YSz
        except:
            continue
        try: 
            for dx in [-1,0,1]:
                for dy in [-1,0,1]:
                    img[x+dx,y+dy,:] = fi
        except: pass

        print x,y

    back = np.median(allspec, 0)

    ff = pf.PrimaryHDU(img.T)
    ff.writeto(outname)

    for cnt, ext in enumerate(exts):
        if ext.xrange is None: continue
        if ext.exptime is None: ext.exptime = 1
        if ext.lamcoeff is None: continue

        ix = np.arange(*ext.xrange)
        l = chebval(ix, ext.lamcoeff)
        s = ext.specw 

        f = interp1d(l, s, fill_value=np.nan, bounds_error=False)
        fi = f(l_grid)/dl_grid 
        fi = fi[0:len(fi):2] + fi[1:len(fi):2] - back


        try:
            x = np.round(Size*np.sqrt(3.) * (ext.Q_ix + ext.R_ix/2)) + XSz
            y = np.round(Size*3./2. * ext.R_ix) + YSz
        except:
            continue
        try: 
            for dx in [-1,0,1]:
                for dy in [-1,0,1]:
                    img[x+dx,y+dy,:] = fi
        except: pass

    ff = pf.PrimaryHDU(img.T)
    ff.writeto("bs_" + outname)


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

        Q_ix: The axial Q coordinate in integral units
        R_ix: The axial R coordinate in integral units

        The relationship of Q/R to X/Y is defined through the pixel mapping
        matrix times the rotation matrix
            P2H = np.array([[np.sqrt(3)/3, -1/3.], [0, 2/3.]]) 
            Rot (22 degree)
    '''

    
    Xs = [None] * len(exts)
    Ys = [None] * len(exts)
    segids = [el.seg_id for el in exts]

    for ext in exts:
        ext.Q_ix = None
        ext.R_ix = None
    
    for ix, ext in enumerate(exts):
        # The X/Y location of a lenslet is based on its
        # trace Y position and where the Halpha wavelength
        # is expected in the X direction.
        if not ext.ok: continue

        Xs[ix] = -999
        Ys[ix] = -999

        if ext.mdn_coeff is not None:
            coeff = ext.mdn_coeff
        elif ext.lamcoeff is not None: 
            coeff = ext.lamcoeff
        else:
            print ext.seg_id, ": ", ext.xrange[0], ext.yrange[0]
            continue

        ixs = np.arange(*ext.xrange)
        LL = chebval(ixs, coeff)
        
        ix_ha = np.nanargmin(np.abs(LL-656.3))

        Xs[ix] = ixs[ix_ha]
        Ys[ix] = np.nanmean(ext.yrange)
        ext.X_pix = Xs[ix]
        ext.Y_pix = Ys[ix]


    Xs = np.array(Xs, dtype=np.float)
    Ys = np.array(Ys, dtype=np.float)
    Xs[Xs != Xs] = -999
    Ys[Ys != Ys] = -999
    dat = np.array([Xs,Ys],dtype=np.float).T

    tree = KDTree(dat)

    ignore, Center = tree.query([1024,1024], 1)
    exts[Center].Q_ix = 0
    exts[Center].R_ix = 0


    def populate_hex(to_populate):
        ''' Breadth-first search 

            For each spaxel in the datacube this piece of code identified
            the relative offset based on the rotation matrix defined
            earlier in the file
        '''
        # Query 7 for the object + six surrounding members
        v = np.array([Xs[to_populate], Ys[to_populate]])
        Dists, Ixs = tree.query(v, 7)
        Tfm = P2H * ROT / np.median(Dists) * 2

        if Dists[0] < 2: 
            Dists = Dists[1:]
            Ixs = Ixs[1:]

        ok = Dists < 70
        Dists = Dists[ok]
        Ixs = Ixs[ok]

        q_this = exts[to_populate].Q_ix
        r_this = exts[to_populate].R_ix

        for nix in Ixs:
            # Search around the current hex via a recrusive call
            # to populate_hex

            nv = np.array([Xs[nix], Ys[nix]])
            D = nv-v
            dp = np.dot(Tfm, D)
            rnd = np.round(dp)

            if rnd[0] == 0 and rnd[1] == 0:
                continue

            if np.abs(rnd[0]) > 1 or np.abs(rnd[1]) > 1:
                continue

            if exts[nix].Q_ix is None:
                exts[nix].Q_ix = q_this + rnd[0]
                exts[nix].R_ix = r_this + rnd[1]
                populate_hex(nix)
            else:
                if (exts[nix].Q_ix != q_this + rnd[0]) or \
                    (exts[nix].R_ix != r_this + rnd[1]):
                    print "collision: "
                    print exts[nix].Q_ix, q_this + rnd[0]
                    print exts[nix].R_ix, r_this + rnd[1]
                    exts[nix].Q_ix = q_this + rnd[0]
                    exts[nix].R_ix = r_this + rnd[1]

    populate_hex(Center)

    # Now convert Q/R to even-Q X/Y
    #

    Qs = np.array([ext.Q_ix for ext in exts], dtype=np.float)
    Rs = np.array([ext.R_ix for ext in exts], dtype=np.float)

    Xs = np.sqrt(3) * (Qs + Rs/2.0)
    Ys = 3/2 * Rs

    # Note 0.633 is plate scale measured on 22 May 2014.
    
    t= np.radians(165.0+45)
    Rot = np.array([[np.cos(t), -np.sin(t)],
                    [np.sin(t),  np.cos(t)]])
    for ix, ext in enumerate(exts):
        p = np.dot(Rot , np.array([Xs[ix], Ys[ix]]))
        ext.X_as = p[0] * -0.633
        ext.Y_as = p[1] * 0.633

    np.save(outname, exts)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Cube.py:

        Convert an extracted file into a data cube

        step is either:
            makecube: To create the data cube (once per night)
            extractcube: To extract the cube (one for each observation)
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('extracted', type=str, help='Extracted file')
    parser.add_argument('--step', type=str, default='make')
    parser.add_argument('--outname', type=str, help='Output cube name')

    args = parser.parse_args()

       
    if args.outname is not None:
        args.outname = args.outname.rstrip('.npy')

    step = args.step
    infile = args.extracted

    if step == 'make':
        print "MAKING"
        ext = np.load(infile)
        cube = extraction_to_cube(ext, outname=args.outname)
    elif step == 'extract':
        print "EXTRACTING"
        ext,meta = np.load(infile)
        QR_to_img(ext, Size=2, outname=args.outname)
    else:
        print "NO STEP TO PERFORM"





