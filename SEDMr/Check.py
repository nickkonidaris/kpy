
import argparse
import numpy as np
import pylab as pl
import pyfits as pf
import sys

 
def checkSpec(specname):

    ss = np.load(specname)[0]

    print "Plotting spectrum in %s" % specname
    print "Extraction radius: %1.2f" % ss['radius_as']
    print ss.keys()

    pl.title("Spectrum in %s" % specname)
    pl.xlabel("Wavelength [nm]")
    pl.ylabel("photon/10 m/nm")
    
    pl.step(ss['nm'], ss['ph_10m_nm'], linewidth=2)
    try: pl.step(ss['skynm'], ss['skyph'])
    except: pl.step(ss['nm'], ss['skyph']*(ss['N_spaxA']+ss['N_spaxB']))

    try: pl.step(ss['nm'], np.sqrt(np.abs(ss['var'])))
    except: pass

    pl.legend(['obj', 'sky', 'err'])
    pl.xlim(360,1100)
    pl.grid(True)
    pl.ioff()
    pl.show()

def checkCube(cubename):
    ''' Plot a datacube for checking'''
    
    cc = np.load(cubename)

    Xs = [c.X_as for c in cc]
    Ys = [c.Y_as for c in cc]
    Ss = [c.trace_sigma for c in cc]

    pl.figure(1)
    pl.scatter(Xs, Ys, marker='H', linewidth=0, s=50, c=Ss, vmin=0.8, vmax= 2)
    pl.title("This should show a regular hexagonal grid of cube positions")
    pl.xlim(-25,25)
    pl.ylim(-25,25)

    pl.colorbar(label='RMS trace width in pixel')
    pl.xlabel("X position [as]")
    pl.ylabel("Y position [as]")
    pl.grid(True)
    pl.ioff()
    pl.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Check.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('--cube', type=str, help='Fine correction path')
    parser.add_argument('--spec', type=str, help='Extracted spectrum file')


    args = parser.parse_args()

    if args.cube is not None:
        checkCube(args.cube)
    if args.spec is not None:
        checkSpec(args.spec)


