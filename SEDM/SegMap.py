'''
SegmentationMap class wraps around several functions in SEDM package and
provides a convenient way to handle SEDM data.
'''

import datetime
import numpy as np
import json
import scipy.io
import matplotlib.pyplot as pl
from matplotlib.backend_bases import KeyEvent 
from matplotlib.backend_bases import PickEvent


import Extract

reload(Extract)

class PositionPicker(object):
    '''Shows the IFU field and allows the user to select an object.

    '''

    index = 0   # index into the file list
    filelist = []
    subtract_sky = True
    picked = []
    outfile = None

    def __init__(self,filelist):
        print "Starting picker GUI"
        self.filelist = filelist
        self.fig = pl.figure(1)
        
        self.fig.canvas.mpl_connect("key_press_event", self)
        self.fig.canvas.mpl_connect("pick_event", self)

        self.index = 0
        self.load()
        self.draw()

    def dump(self):
        """Write status to file"""

        dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print "Writing to: %s" % self.outfile
        str = json.dumps({"outfile": self.outfile,
                        "infile": self.filelist[self.index],
                        "spec": self.picked,
                        "status": self.status,
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
        self.SM = SegmentationMap(self.filelist[self.index])

        self.outfile = self.filelist[self.index].rstrip("fits_SI.mat") + \
                        ".coords.json"
        print "Loaded"

    def __call__(self, event):
        '''Event call handler for Picker gui.'''

        if event.name == 'pick_event':
            self.picked = self.SM.OK[event.ind[0]]
            pl.title("[1] Clicked on %s. [2] Press q to exit." % event.ind)

        elif event.name == 'key_press_event':
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
                print "Substract sky: %s" % self.subtract_sky
                self.subtract_sky = not self.subtract_sky
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

    
        x,y,v = Extract.segmap_to_img(self.SM.SegMap, sky_spec=sky_spec)
        self.Xs = x
        self.Ys = y
        self.Values = v


        pl.ion()
        pl.figure(1, figsize=(9,8))
        pl.clf()
        pl.xlim(-12,17)
        pl.ylim(-7,27)
        pl.title("[1] Click on an object. [2] Label as [g]ood, [u]nsure, [n]othing visible. [3] Goto [n]ext/[p]rev object")

        pl.scatter(self.Xs,
                    self.Ys,
                    c=self.Values,
                    s=60,
                    picker=0.5,
                    marker='h')

        
        pl.colorbar()

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

    def __init__(self, segmap=None):
        ''' Loads the segmap, creates the KD Tree.

        Args:
            segmap: Either a filename or a segmentation map array '''
        
        if type(segmap) == str:
            print "Loading: %s" % segmap
            mat = scipy.io.loadmat(segmap)
            self.SegMap = mat['SegmentsInfo']
        else:
            self.SegMap = segmap

        self.KT, self.OK = Extract.segmap_to_kdtree(self.SegMap)
    
    def subtract(self, B):
        '''Subtracts B.SegMap from current SegMap.

        sm.substract(B) has side effects and overwrites the segmap in self

        Args:
            B: The SegmentationMap object to subtract off

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

            d = self.SegMap[i]['SpexSpecFit'][0][0] \
                    -  B.SegMap[i]['SpexSpecFit'][0][0]

            self.SegMap[i]['SpexSpecFit'][0][0] = d

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
            sky_spec=sky_spec)

        fig = pl.figure(figure_number, figsize=(9,8))
        pl.xlim(-12, 17)
        pl.ylim(-7, 27)
        pl.scatter(x,y,c=v, s=60, picker=0.5, marker='h')

    def select_circle(self, center_xs, distance=2):
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

        if center_xs is None:
            return None

        ix = self.OK[center_xs[0]]
        X = self.SegMap[ix]['OnSkyX'][0][0][0]
        Y = self.SegMap[ix]['OnSkyY'][0][0][0]

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

       
        res = Extract.sky_median(self.KT, self.SegMap, ixmap=self.OK,
            pixel_shift=self.pixel_shift)

        return res

