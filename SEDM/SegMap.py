'''
SegmentationMap class wraps around several functions in SEDM package and
provides a convenient way to handle SEDM data.
'''

import numpy
import scipy.io
import pylab as pl
import matplotlib.pyplot as plt


import Extract

reload(Extract)

picked_points = None
artist = None
    
def pick_handler(event):
    global picked_points,artist

    mouseevent = event.mouseevent
    artist = event.artist

    picked_points = event.ind
    print "Picked: %s" % picked_points

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
            mat = scipy.io.loadmat(segmap)
            self.SegMap = mat['SegmentsInfo']
        else:
            self.SegMap = segmap

        self.KT, self.OK = Extract.segmap_to_kdtree(self.SegMap)
    
    
    def draw(self, figure_number=1, minl=450, maxl=800, 
        subtract_sky=False):
        '''Draw a figure showing the segmentation map cube
        
        Args:
            figure_number: Plot figure number
            minl/maxl: Minimum/Maximum wavelength to sum over
        
        Returns:
            Nothing
            
        Side Effects:
            Plots a new figure(figure_number) with the data cube.
            global picked_points will be updated with the picked points'''

        global picked_points
        picked_points = None

        sky_spec = None
        if subtract_sky:
            if self.sky_spectrum[0] is None:
                raise Exception("Request to subtract sky; however, no sky has been measured")
            sky_spec = self.sky_spectrum

        x,y,v = Extract.segmap_to_img(self.SegMap, minl=minl, maxl=maxl, 
            sky_spec=sky_spec)

        fig = plt.figure(figure_number, figsize=(8,7))
        ax = fig.add_subplot(111)
        plt.clf()
        plt.scatter(x,y,c=v,s=60, picker=True, marker='h')
        plt.ylim(-7,27)
        plt.xlim(-12,17)
        fig.canvas.mpl_connect('pick_event', pick_handler)
        plt.show()

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
                self.spectrum_near_position(X,Y, distance, 
                    pixel_shift=self.pixel_shift),
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

