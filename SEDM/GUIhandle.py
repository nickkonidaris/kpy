'''
SegmentationMap class wraps around several functions in SEDM package and
provides a convenient way to handle SEDM data.
'''

import datetime
import numpy as np
import json
import pyfits
import scipy.io
import matplotlib.pyplot as pl
from matplotlib.backend_bases import KeyEvent 
from matplotlib.backend_bases import PickEvent

import Disp


import Extract

reload(Extract)
reload(Disp)

class PositionPicker(object):
    '''Shows the IFU field and allows the user to select an object.

    '''

    index = 0   # index into the file list
    filelist = []
    subtract_sky = True
    picked = []
    outfile = None
    olines  = [372.7, 486.1, 500.7, 656.3]
    tlines = [761.5, 589.0, 557.7, 435.8, 519.9, 630.0]
    pixel_shift = 0
    positions=None
    qecurve = None
    

    def __init__(self,filelist,positions=('OnSkyX','OnSkyY'),
        qefunction=None):
        print "Starting picker GUI"
        self.filelist = filelist
        self.positions=positions
        self.fig = pl.figure(1)

        self.fig.canvas.mpl_connect("key_press_event", self)
        self.fig.canvas.mpl_connect("pick_event", self)
        self.fig.canvas.mpl_connect("button_press_event", self)

        self.fig2 = pl.figure(2)
        self.fig2.canvas.mpl_connect("key_press_event", self)
        self.fig2.canvas.mpl_connect("pick_event", self)

        self.qecurve = qefunction
        self.index = 0
        self.load()
        self.draw()

    def dump(self):
        """Write status to file"""

        dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        X,Y = self.picked
        print "Writing to: %s" % self.outfile
        str = json.dumps({"outfile": self.outfile,
                        "infile": self.filelist[self.index],
                        "spec": (X,Y),
                        "status": self.status,
                        "pixel_shift": self.pixel_shift,
                        "when": dt}, indent=4)

        try:
            f=open(self.outfile, "w")
            f.write(str)
            f.close()
        except Exception as e:
            raise("Could not write %s: %s" % (self.outfile, e))

    def load(self):
        """Load the Segmentation map"""

        pl.figure(1)
        pl.clf()
        if self.index < 0: self.index =0
        if self.index > len(self.filelist): index = len(self.filelist)-1

        self.picked = []
        self.pixel_shift = 0
        self.SM = SegmentationMap(self.filelist[self.index],
            positions=self.positions)

        try: reg_txt = Disp.ds92(self.SM.SegMap)
        except: pass

        try:
            outreg = self.filelist[self.index].rstrip("fits_SI.mat") + \
                        ".reg"
            f = open(outreg, "w")
            f.write(reg_txt)
            f.close()
        except:
            print "Could not write reg file"

        self.outfile = self.filelist[self.index].rstrip("fits_SI.mat") + \
                        ".coords.json"
        print "Loaded"

    def draw_spectrum(self):
        ''' Draw the spectrum in figure(2) '''


        pl.figure(2)
        pl.clf()
        if self.picked == []: return
        self.SM.pixel_shift = self.pixel_shift

        X,Y = self.picked
        print "Extracting at: %s,%s" % (X,Y)
        sky_spec = self.SM.spectrum_in_annulus(X,Y)
        wave = sky_spec["wave_nm"]
        sky = (wave, sky_spec["spec_adu"]/sky_spec["num_spec"])

        obj_spec = self.SM.spectrum_near_position(X,Y, onto=wave)

        pl.xlabel("wavelength [nm]")
        spec = obj_spec["spec_adu"]

        if self.subtract_sky:
            spec -= sky[1]*obj_spec["num_spec"]

        self.sky_spec = sky
        self.obj_spec = (wave, spec)

        if self.qecurve is None: correction = 1.0
        else:
            correction = 1/self.qecurve(wave*10)
            correction[correction<.01] = 1.0
            correction[correction>100]=1.0

        pl.step(wave, spec*correction)
        pl.step(wave, sky[1]*obj_spec["num_spec"]*correction,'r')

        pl.legend(["object","sky"])
        pl.xlim(350,1000)

        for line in self.olines:
            pl.axvline(line)
        for line in self.tlines:
            pl.axvline(line, color='r')

    def handle_shift(self, xdata, ydata):

        if (xdata < 360) or (xdata > 1000): return

        lines = np.concatenate((self.olines, self.tlines))

        print lines
        delts = (lines - xdata)
        ix = np.argmin(np.abs(delts))
        print "Closest to %f" % lines[ix]
        
        line = lines[ix]
        delt = delts[ix]
        wave = self.sky_spec[0]

        wix = np.nanargmin(np.abs(wave-line))
        dw = wave[wix]-wave[wix-1]
        print "Delt: {0}, dw: {1}".format(delt, dw)
        self.pixel_shift += delt/dw

        print "pixel shift is: {0}".format(self.pixel_shift)
        self.draw_spectrum()

    def __call__(self, event):
        '''Event call handler for Picker gui.'''

        if event.name == 'button_press_event':
            self.picked = (event.xdata, event.ydata)
            self.draw_spectrum()
            
        elif event.name == 'key_press_event':
            if event.key == '\\':
                print "Shifting"
                self.handle_shift(event.xdata, event.ydata)

            if event.key == 'n': 
                print "next"
                self.index += 1
                self.load()
                self.draw()
            if event.key == 'p': 
                print "prev"
                self.index -= 1
                self.load()
                self.draw()
            if event.key == '-':
                self.subtract_sky = not self.subtract_sky
                print "Substract sky: %s" % self.subtract_sky
                self.draw()
            if event.key == "u":
                self.status = "unsure"
                self.dump()
            if event.key == "b":
                self.status = "bad"
                self.dump()
            if event.key == "o":
                self.status = "ok"
                self.dump()

            if event.key == 'h':
                print """Help---
n - next
p - prev
- - subtract sky
u - unsure: there are targets visible, not sure which is correct.
b - bad: nothing visible
o - ok: target visible
"""
            print event.key
            

    def draw(self):
        if self.subtract_sky:
            sky = self.SM.sky_median()
            sky_spec = sky['wave_nm'], sky['spec_adu']/sky['num_spec']
        else:
            sky_spec = None

    
        x,y,v = Extract.segmap_to_img(self.SM.SegMap, sky_spec=sky_spec,
            minl=500, maxl=700,positions=self.positions)
        self.Xs = x
        self.Ys = y
        self.Values = v

        pl.figure(1, figsize=(9,8))

        pl.ion()
        pl.clf()
    
        try: fname = self.filelist[self.index].split("/")[-1]
        except: fname = "???"
        try: 
            fitsname = self.filelist[self.index]
            fitsname = fitsname.replace("shrunk_","").rstrip("_SI.mat")

            header = pyfits.getheader(fitsname)
            name = header['OBJECT']
        except: name = "???"

        if self.positions[0] == 'OnSkyX': diam = 0.5
        else: diam = 1


        pl.title("{0}/{1}".format(name, fname))
        pl.scatter(self.Xs,
                    self.Ys,
                    c=self.Values,
                    s=40,
                    picker=diam,
                    marker='h')
        pl.colorbar()


        self.draw_spectrum()
        

