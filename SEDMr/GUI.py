
import NPK.Atmosphere as Atm
import NPK.PlotHelp as PH
import datetime
import numpy as np
import os
import pyfits
import scipy.io
import matplotlib.pyplot as pl
from matplotlib import gridspec
from matplotlib.backend_bases import KeyEvent 
from matplotlib.backend_bases import PickEvent
from scipy.interpolate import interp1d
from matplotlib.widgets import Cursor

import scipy, scipy.spatial

from numpy.polynomial.chebyshev import chebfit, chebval


import SEDMr.Spectra as SS


class MouseCross(object):
    ''' Draw a cursor with the mouse cursor '''

    def __init__(self, ax, radius_as=3, **kwargs):
        self.ax = ax

        radius_pix = np.abs((ax.transData.transform((radius_as, 0)) -
            ax.transData.transform((0,0)))[0])
        print "%s arcsec is %s pix" % (radius_as, radius_pix)
        self.line, = self.ax.plot([0], [0], visible=False, 
            marker=r'$\bigodot$', markersize=radius_pix*2, color='red', **kwargs)

    def show_cross(self, event):
        if event.inaxes == self.ax:
            self.line.set_data([event.xdata], [event.ydata])
            self.line.set_visible(True)
        else:
            self.line.set_visible(False)


        pl.draw()

class PositionPicker(object):
    ''' This class is used to select an extraction point in a data cube '''

    spectra = None
    Xs = None
    Ys = None
    Vs = None
    pointsize = None
    picked = None
    radius_as = None

    def __init__(self, spectra=None, figure=None, pointsize=55, bgd_sub=False, radius_as=3):
        ''' Create spectum picking gui.

        Args:
            spectra: SEDMr.Spectra object
            figure: Figure to draw to [default is None or open a new figure
                window '''

        self.spectra = spectra
        self.pointsize = pointsize

        self.Xs, self.Ys, self.Vs = spectra.to_xyv()

        if bgd_sub:
            self.Vs -= np.median(self.Vs)

        pl.ioff()
        self.figure = pl.figure(1)

        self.radius_as = radius_as

        self.figure.canvas.mpl_connect("button_press_event", self)
        self.draw_cube()

    def draw_cube(self):
        pl.scatter(self.Xs, self.Ys, c=self.Vs, s=self.pointsize, linewidth=0)
            
        pl.ylim(-20,20)
        pl.xlim(-20,20)
        pl.colorbar()

        c = Cursor(self.figure.gca(), useblit=True)

        cross = MouseCross(self.figure.gca(), radius_as=self.radius_as)
        self.figure.canvas.mpl_connect('motion_notify_event', 
            cross.show_cross)
        pl.show()

    def __call__(self, event):
        '''Event call handler for Picker gui.'''
        
        if event.name == 'button_press_event':
            print event.xdata, event.ydata
            self.picked = (event.xdata, event.ydata)
            pl.close(self.figure)
            


class WaveFixer(object):
    ''' This class is used to fix bad wavelength solutions '''

    cube = None # Raw data cube spectra
    KT = None # KDTree object
    X1 = []
    X2 = []
    Y1 = []

    pointsize = None
    picked = None

    state = "Display"

    fig = None
    ax_cube  = None
    ax_spec  = None


    def __init__(self, cube=None, figure=None, pointsize=55, bgd_sub=False, radius_as=3):
        ''' Create spectum picking gui.

        Args:
            cube: Data cube list
            figure: Figure to draw to [default is None or open a new figure
                window '''

        self.actions = {"m": self.mode_switch}

        self.cube = cube
        self.pointsize = pointsize

        for ix, s in enumerate(self.cube):
            if s.xrange is None: continue
            xs = np.arange(s.xrange[0], s.xrange[1], .1)

            if s.lamcoeff is not None:
                lls = chebval(xs, s.lamcoeff)
                ha1 = np.argmin(np.abs(lls - 656.3))/10.0 + s.xrange[0]

            if s.mdn_coeff is not None:
                lls = chebval(xs, s.mdn_coeff)
                ha2 = np.argmin(np.abs(lls - 656.3))/10.0 + s.xrange[0]

            self.X1.append(ha1)
            self.X2.append(ha2)
            self.Y1.append(s.yrange[0])


        self.X1, self.X2, self.Y1 = map(np.array, [self.X1, self.X2,
            self.Y1])

        OK = (np.abs(self.X1 - self.X2) < 2) & np.isfinite(self.X2) & \
            np.isfinite(self.Y1)

        self.good_cube = self.cube[OK]
        locs = np.array([self.X2[OK], self.Y1[OK]]).T
        self.KT = scipy.spatial.KDTree(locs)

        assert(len(locs) == len(self.good_cube))

        # Setup drawing
        gs = gridspec.GridSpec(1,2, width_ratios=[1,2.5])
        self.fig = pl.figure(1, figsize=(22,5.5))
        self.ax_cube = pl.subplot(gs[0])
        self.ax_spec = pl.subplot(gs[1])

        self.ax_cube.set_xlim(-100, 2200)
        self.ax_cube.set_ylim(-100, 2200)

        pl.ion()
        self.draw_cube()

        print "Registering events"
        self.fig.canvas.mpl_connect("button_press_event", self)
        self.fig.canvas.mpl_connect("key_press_event", self)

        pl.show()


    def mode_switch(self):
        ''' Toggle operating mode between Display and Select '''

        if self.state == "Display": self.state = "Select"
        else: self.state = "Display"

        print self.state
        self.draw_cube()

    def draw_spectra(self):
        ''' Draw nearest spectra '''
        if self.picked is None: return

        print "Drawing spectra"
        #xl = self.ax_spec.get_xlim()
        #yl = self.ax_spec.get_ylim()
        self.ax_spec.cla()
        #self.ax_spec.set_xlim(xl)
        #self.ax_spec.set_ylim(yl)

        
        x,y = self.X2[self.picked], self.Y1[self.picked]
        objs = self.KT.query_ball_point((x,y), 70)

        print "Query around %s found: %s" % ((x,y), objs)

        spec = self.cube[self.picked]
        ix = np.arange(*spec.xrange)
        fiducial_ll = chebval(ix, spec.lamcoeff)
        #self.ax_spec.plot(ix-spec.xrange[0], fiducial_ll, linewidth=3)
        self.ax_spec.step(fiducial_ll, spec.spec, linewidth=3)

        for spec in self.good_cube[objs]:
            try: 
                ix = np.arange(*spec.xrange)
                ll = chebval(ix, spec.lamcoeff)
                #self.ax_spec.plot(ix-spec.xrange[0], ll-fiducial_ll)
                self.ax_spec.step(ll, spec.spec)
            except: pass


        #self.ax_spec.set_ylim(-30,30)
        self.ax_spec.set_xlim(370,700)
        self.fig.show()


    def draw_cube(self):
        ''' Draw the data cube '''
        print "drawing cube"
        
        # Draw cube

        xl = self.ax_cube.get_xlim()
        yl = self.ax_cube.get_ylim()
        self.ax_cube.cla()
        self.ax_cube.set_xlim(xl)
        self.ax_cube.set_ylim(yl)

        self.ax_cube.plot(self.X2, self.Y1, 'o',
            markersize=8, marker='h', linewidth=0)

        bad = np.abs(self.X1-self.X2) > 4
        self.ax_cube.plot(self.X2[bad], self.Y1[bad], 'ro',
            markersize=7, marker = 'h', linewidth=0)

        if self.picked is not None:
            print self.X2[self.picked], self.Y1[self.picked]

            self.ax_cube.plot([self.X2[self.picked]], 
                [self.Y1[self.picked]], 'o', ms=12, color='yellow',
                alpha=0.4, visible=True)

        tit = "State: %s | Press ? for help" % self.state
        self.ax_cube.set_title(tit)

        self.fig.show()


    def handle_button_press(self, event):
        

        if event.inaxes == self.ax_cube:
            '''Clicked In Data Cube Display'''

            if self.state == 'Display':
                ''' Display state (not pick state, ignore) '''
                return

            dists = np.abs(self.X2 - event.xdata) + \
                np.abs(self.Y1 - event.ydata)

            ix = np.nanargmin(dists)
            print dists[ix]

            if dists[ix] > 20: self.picked = None
            else: self.picked = ix

            self.draw_cube()
            self.draw_spectra()





    def __call__(self, event):
        '''Event call handler for Picker gui.'''
        
        print (event.name)

        if event.name == 'pick_event':
            import pdb
            pdb.set_trace()

        elif event.name == 'button_press_event':
            ''' Note order of if statement to skip button over pick event'''

            self.handle_button_press(event)
            
        elif event.name == 'key_press_event':
            key = event.key

            if key in self.actions:
                to_call = self.actions[key]
                to_call()

            if key == "?":
                for k,v in self.actions.iteritems():
                    print "%s: %s" % (k,v.__doc__)





