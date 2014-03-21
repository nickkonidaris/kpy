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
        

class SegmentationMap(object):
    KT = None # KD Tree
    OK = None # list of segmap elements with np.isfinite OnSkyX and OnSkyY
        # positions. Note that the KT index values return an index into
        # OK and not SegMap. Any KT Queries will require something like:
        # ixs = KT.query_(x,y ..)
        # SegMap[OK[ixs[i]]]
    SegMap = None # The segmentation map
    sky_spectrum = (None, None) # Sky spectrum (Wavelength, Spectrum)
        # sky spectrum is a per-spaxel spectrum
    pixel_shift = 0 # The number of pixels to translate the spectra. pixel_shift
        # is used for flexure correction and is measured from the spectra
        # themselves.

    positions=None # Use OnSky? or Mean? from the SegmentationMap
    signal_field=None # Use SpexSpecFit or SpexSpecCenter

    def __init__(self, segmap=None, positions=('OnSkyX', 'OnSkyY'),
        signal_field='SpexSpecFit'):
        ''' Loads the segmap, creates the KD Tree.

        Args:
            segmap: Either a filename or a segmentation map array '''
        
        if type(segmap) == str:
            print "Loading: %s" % segmap
            mat = scipy.io.loadmat(segmap)
            self.SegMap = mat['SegmentsInfo']
        else:
            self.SegMap = segmap

        self.positions=positions
        self.signal_field=signal_field
        self.KT, self.OK = Extract.segmap_to_kdtree(self.SegMap,
            positions=positions,signal_field=signal_field)
    
    def subtract(self, B):
        '''Subtracts B.SegMap from current SegMap.

        sm.substract(B) has side effects and overwrites the segmap in self

        Args:
            B: The SegmentationMap object to subtract off
            signal_field: String containing the field to extract signal from
                defaults to SpexSpecFit, could also be SpexSpecCenter
            

        Returns:
            Nothing, but has side effects.'''

    
        if len(self.SegMap) != len(B.SegMap):
            raise Exception("Mismatched segmentation maps")

        for i in xrange(len(self.SegMap)):

            try:
                w = self.SegMap[i]['WaveCalib'][0][0]
            except:
                continue

            ok = np.isfinite(w) & np.isfinite(B.SegMap[i]['WaveCalib'][0][0])
            if np.any(w[ok] != B.SegMap[i]['WaveCalib'][0][0][ok]):
                raise Exception("Mismatched segmentation maps, wave")

            d = self.SegMap[i][self.signal_field][0][0] \
                    -  B.SegMap[i][self.signal_field][0][0]

            self.SegMap[i][signal_field][0][0] = d

    def draw(self, figure_number=1, minl=450, maxl=800, 
        subtract_sky=False, outfile=None):
        '''Draw a figure showing the segmentation map cube
        
        Args:
            figure_number: Plot figure number
            minl/maxl: Minimum/Maximum wavelength to sum over
            outfile: The path to the output file
        
        Returns:
            Nothing
            
        Side Effects:
            Plots a new figure(figure_number) with the data cube.
        '''

        sky_spec = None
        if subtract_sky:
            if self.sky_spectrum[0] is None:
                raise Exception("Request to subtract sky; however, no sky has been measured")
            sky_spec = self.sky_spectrum

        x,y,v = Extract.segmap_to_img(self.SegMap, minl=minl, maxl=maxl, 
            sky_spec=sky_spec,positions=self.positions, 
            signal_field=self.signal_field)

        fig = pl.figure(figure_number, figsize=(9,8))
        if self.positions[0]=='MeanX':
            pl.xlim(-100,2200)
            pl.ylim(-100,2200)
        else:
            pl.xlim(-12, 17)
            pl.ylim(-7, 27)
        pl.scatter(x,y,c=v, s=60, picker=0.5, marker='h')

    def select_circle(self, center_xs, distance=2, 
        positions=('OnSkyX','OnSkyY')):
        '''
        Selects and reutrns the object spectra near spectrum index `center_xs`

        Args:
            [center_xs]: center_xs[0] is the spectrum index to search around.
                all other values in center_xs ignored.
            distance: Radial distance to search around in arcsecond.
        Returns:
            List of spectra:
                [[Wavelength [nm], Counts [ADU], spectrum index]]

            Note the spectra are on a different wavelength grid.
        '''

        if self.positions[0]=='MeanX':
            if distance<60: distance*=60

        if center_xs is None:
            return None

        ix = self.OK[center_xs[0]]
        X = self.SegMap[ix][positions[0]][0][0][0]
        Y = self.SegMap[ix][positions[1]][0][0][0]

        print X,Y

        spectra = self.spectra_near_position( X,Y,
                                        distance=distance)

        return spectra

    def spectra_near_position(self, X, Y, distance=2):
        '''Returns all spectra within distance of (X,Y)

        Notice spectra_* methods return all spectra, spectrum_* methods
            convert spectra into a single spectrum.

        Args:
            X,Y: The X/Y position of the spectrum to extract in as
            distance: The extraction radius in arcsecond

        Returns a list of (wavelength, spectra, index). E.g.:
            [[array(365...1000) , array(50 .. 63), 950] ..]

        Example:
            import SegMap as SM
            sm = SM.SegmentationMap("/path/b_ifu20130808_23_08_44.fits_SI.mat")
            spec = sm.spectra_near_position(5, 5)
            print len (spec)
                >> 34
            ll, ss = spec[0][0], spec[0][1]
            plot(ll, ss) # Plots the spectrum
        '''

        if self.positions[0]=='MeanX':
            if distance<60: distance*=60

        return Extract.spectra_near_position(self.KT, self.SegMap,
                        X,Y, distance=distance, ixmap=self.OK, 
                        pixel_shift=self.pixel_shift)

    def spectrum_near_position(self, X, Y, distance=2, onto=None, 
                                sky_spec=None):
        '''Interpolate spectra in a circle radius distance arround X,Y

        See spectrum_near_position

        Args: 
            X,Y: The X/Y position of the central spectrum
            small, large: the small and large radius of extraciton
            onto: The wavelength grid to interpolate onto, None to ignore.
            sky_spec: The sky spectrum to subtract

        Returns:
            {'wave_nm': wavelength of sky spectrum, 
                'spec_adu', the median sky spectrum, 
                'all_spec': and a matrix of spectra,
                'num_spec': The number of spectra}'''

        if self.positions[0]=='MeanX':
            if distance<60: distance*=60

        return Extract.interp_and_sum_spectra(
                self.spectra_near_position(X,Y, distance),
                onto=onto, sky_spec=sky_spec)
        
    def spectra_in_annulus(self, X, Y, small=4, large=6):
        '''Returns all spectra in the annulus between radius small and large.

        Notice spectra_* methods return all spectra, spectrum_* methods
            convert spectra into a single spectrum.

        Args:
            X,Y: The X/Y position of the central spectrum
            small,large: The small and large radius of extraction

        Returns:
            List of spectra, see spectra_near_position'''

        if self.positions[0]=='MeanX':
            if small<60: small *= 60
            if large<60: large *= 60
            
        return Extract.spectra_in_annulus(self.KT, self.SegMap,
                X, Y, small=small, large=large, ixmap=self.OK, 
                pixel_shift=self.pixel_shift)

    def spectrum_in_annulus(self, X, Y, small=4, large=6, onto=None):
        '''Returns the interpolated spectrum in an annulus arround X,Y

        See spectra_in_annulus

        Args: 
            X,Y: The X/Y position of the central spectrum
            small, large: the small and large radius of extraciton
            onto: The wavelength grid to interpolate onto, None to ignore.

        Returns:
            {'wave_nm': wavelength of sky spectrum, 
                'spec_adu', the median sky spectrum, 
                'all_spec': and a matrix of spectra,
                'num_spec': The number of spectra}'''

        if self.positions[0]=='MeanX':
            if small<60: small *= 60
            if large<60: large *= 60

        return Extract.interp_and_sum_spectra(
                self.spectra_in_annulus(X,Y, small, large))
                    
                
            

    def sky_median(self, X=2, Y=10, distance=5):
        '''Generates a sky spectrum from the median of a large sky area.

        Args:
            X/Y: the X/Y position of the center coordinate
            distance: the radius around the center coordinate

        Returns:
            {'wave_nm': wavelength of sky spectrum, 
                'spec_adu', the median sky spectrum, 
                'all_spec': and a matrix of spectra,
                'num_spec': The number of spectra}

        Side Effects:
            Stores the returned values into the class'''

       
        if self.positions[0]=='MeanX':

            if X<60:
                X = 1024
                Y=1024
            if distance<60: distance*= 60

        res = Extract.sky_median(self.KT, self.SegMap, ixmap=self.OK,
            pixel_shift=self.pixel_shift,X=X,Y=Y,distance=distance)

        return res

