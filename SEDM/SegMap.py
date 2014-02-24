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
        spectra = Extract.spectra_near_position(self.KT, self.SegMap,
                                        X,Y,
                                        distance=distance)

        return spectra

    def spectra_near_position(
        # NEXT STEP IS HERE

    def sky_median(self, X=2, Y=10, distance=7):
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

       
        res = Extract.sky_median(self.KT, self.SegMap)

        return res
            
        




