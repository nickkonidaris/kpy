import argparse 
import pylab as pl
import numpy as np
import sys
import os
from scipy.interpolate import interp1d
from scipy.interpolate import splrep,splev
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.ndimage.filters import gaussian_filter

import NPK.Standards as Stds

def near_balmer(lam, Rwindow=.02):
    lines = [656.3, 486.1, 434.1, 410.2, 397.0]

    if np.min(np.abs(lam-lines)/lines) < Rwindow:
        return True
    else:
        return False

class normalizer:
    """
    An interactive spline fitting routine for determining the continuum
     level of a spectrum.

    """

    def __init__(self, wave, flux, sss, window=15.0, name=None, filename=None):
        self.wave = wave
        self.flux = flux
        self.continuum = None
        self.cfunc = None
        self.locations = []
        self.values = []
        self.std_star_fun = sss
        fig = pl.figure()
        self.ax = pl.gca()
        # the width of the median window taken at each point
        self.winwidth = window/2.0

        self.ax.step(self.wave,self.flux,where='mid',label='spectrum')
        pl.grid(True)
        if name != None:
            pl.title(name)
        self.filename = filename

        # Connect the different functions to the different events
        fig.canvas.mpl_connect('key_press_event',self.ontype)
        fig.canvas.mpl_connect('button_press_event',self.onclick)
        fig.canvas.mpl_connect('pick_event',self.onpick)

        print '*'*20
        print ' L click to define spline points'
        print ' R click on a point to remove it'
        print ' press enter at any point to refit continuum'
        print ' after continuum is fit, press "n" to normalize'
        print ' press "r" to reset'
        print ' press "w" to write normalized spectrum out to file'
        print 'When finished, the continuum is accessible as "<normalizer>.continuum"'
        print 'The continuum function is accessible as "<normalizer>.cfunc"'
        print '*'*20

        # Now put some control points down
            
        for lam in np.logspace(np.log10(350), np.log10(1100), 20):
            if near_balmer(lam): continue

            window = ((lam-self.winwidth)<=self.wave) &\
                      (self.wave<=(lam+self.winwidth))
            y = np.median(self.flux[window])

            if y != y: 
                continue

            self.ax.plot(lam,y,'rs',ms=10,picker=5,label='cont_pnt')

        pl.show() # show the window


    def onclick(self, event):
        # when none of the toolbar buttons is activated and the user clicks in the
        # plot somewhere, compute the median value of the spectrum in a 10angstrom
        # window around the x-coordinate of the clicked point. The y coordinate
        # of the clicked point is not important. Make sure the continuum points
        # `feel` it when it gets clicked, set the `feel-radius` (picker) to 5 points
        toolbar = pl.get_current_fig_manager().toolbar
        if event.button==1 and toolbar.mode=='':
            window = ((event.xdata-self.winwidth)<=self.wave) &\
                      (self.wave<=(event.xdata+self.winwidth))
            y = np.median(self.flux[window])
            if y != y: 
                print "clicked out of range"
                return

            if near_balmer(event.xdata):
                self.ax.plot(event.xdata,y,'bs',ms=10,picker=5,label='cont_pnt')
                print "Blue point near balmer line"
            else:
                self.ax.plot(event.xdata,y,'rs',ms=10,picker=5,label='cont_pnt')

        pl.draw()

    def onpick(self, event):
        # when the user right clicks on a continuum point, remove it
        if event.mouseevent.button==3:
            if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
                event.artist.remove()

    def ontype(self, event):
        # when the user hits enter:
        # 1. Cycle through the artists in the current axes. If it is a continuum
        #    point, remember its coordinates. If it is the fitted continuum from the
        #    previous step, remove it
        # 2. sort the continuum-point-array according to the x-values
        # 3. fit a spline and evaluate it in the wavelength points
        # 4. plot the continuum
        if event.key=='enter':
            cont_pnt_coord = []
            for artist in self.ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                    cont_pnt_coord.append(artist.get_data())
                elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                    artist.remove()

            cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
            sort_array = np.argsort(cont_pnt_coord[:,0])
            x,y = cont_pnt_coord[sort_array].T
            self.locations = x
            self.values = y
            spline = splrep(x,y,k=3)
            self.continuum = splev(self.wave,spline)
            self.cfunc = lambda w: splev(w, spline)
            pl.plot(self.wave,self.continuum,'r-',lw=2,label='continuum')

        # when the user hits 'n' and a spline-continuum is fitted, normalise the
        # spectrum
        elif event.key=='n':
            if self.continuum is not None:
                self.ax.cla()
                std_spec = self.std_star_fun(self.wave)

                # Convole spectrum to SEDM Resolution
                std_spec = gaussian_filter(std_spec, 1.5)
                self.correction = std_spec / self.continuum

                self.ax.plot(self.wave,self.flux*self.correction,'k-',label='normalised')

        # when the user hits 'r': clear the axes and plot the original spectrum
        elif event.key=='r':
            self.continuum = None
            self.locations = []
            self.values = []
            self.ax.cla()
            self.ax.plot(self.wave,self.flux,'k-')

        elif event.key == 'q':
            sys.exit()
        # when the user hits 'w': if the normalised spectrum exists, write it to a
        # file.

        elif event.key=='w':
            for artist in self.ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='normalised':
                    data = np.array(artist.get_data())
                    if self.filename == None:
                        self.filename = 'normalised_spec.flm'
                    
                    results = [{
                        'filename': self.filename,
                        'locations_nm': self.locations,
                        'values': self.values,
                        'continuum': self.continuum,
                        'correction': self.correction,
                        'wave': self.wave
                    }]
                    np.save(self.filename, results)
                    print('Saved to file:',self.filename)
                    break
        pl.draw()


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description=\
    '''SpecNorm.py measures the continuum of a standard star spectrum''')

    parser.add_argument('spectrum', type=str, help='Numpy file of spectrum extraction')
    parser.add_argument('--window', type=float, help='Width (nm) to median the splien points over')
    parser.add_argument('--spline_points', type=str, help='Path to ascii text file specifying where to put spline knots')
    parser.add_argument('--std', type=str, help='Name of standard star (see NPK.Standard for acceptable names)')
    parser.add_argument('--outname', type=str, help='Name of output file')

    args = parser.parse_args()

    
    sp = np.load(args.spectrum)[0]
    w = sp['nm']
    s = sp['ph_10m_nm']

    if args.std is None: 
        sss = lambda x: 1
    elif args.std not in Stds.Standards:
        print "SpecNorm.py: Could not find the standard named %s" % args.std
        sss = lambda x: 1
    else:
        lstd = Stds.Standards[args.std][:,0]
        sstd = Stds.Standards[args.std][:,1]
        lstd /= 10.0
        fun = interp1d(lstd, sstd, bounds_error=False, fill_value=np.nan)


    if args.outname is None:
        path = os.path.dirname(args.spectrum)
        filename = os.path.basename(args.spectrum)
        fname = os.path.join(path, "corr_"+filename)
    else:
        fname = args.outname

    nr = normalizer(w,s, fun, filename=fname)

