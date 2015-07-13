import NPK.Atmosphere as Atm
import NPK.PlotHelp as PH
import datetime
import numpy as np
import os
import pyfits
import scipy.io
import matplotlib.pyplot as pl
from matplotlib.backend_bases import KeyEvent 
from matplotlib.backend_bases import PickEvent
from scipy.interpolate import interp1d
from matplotlib.widgets import Cursor

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
        self.figure = pl.figure(figure)

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
            

