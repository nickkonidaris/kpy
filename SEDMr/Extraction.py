'''

SED Machine Extraction class
(c) 2014 Nick Konidaris


This class is an abstract wrapper around the extraction

'''


class Extraction():
    '''SED machine spectral extraction class.

    Generally there are several thousand such Extractions that are held for
    a single 2k x 2k frame. 

    Args:
        seg_id (int): Number from 1 to the maximum segment number. The segment
            number is generated by sextractor.
        ok (bool): Flag indicating the success of wavelength solutions
        xrange([int,int]): From sextractor, the beginning and end of the segment.
            Usually ~ 250 pixels long.
        yrange([int,int]): From sextractor, the top and bottom of the segment.
            Usually less than ~4 pixels.
        poly ([float]): Coefficients to the polynomial function tracing the
            ridge of each spectrum. Is converted to a polynomial function via
            np.poly1d. The function parameter is a pixel integer starting
            from xrange[0]
        spec ([float]): The extracted 1d spectrum.
        hg_lines ({wavelength: pixel}): An association of mercury lamp wavelength
            and pixel position. These values come from sextractor
        hgcoef ([float]): Chebyshev polynomial coefficients to the best-fit
            mercury lamp spectrum. Used as an intermediate step. Do not use
            in general.
        lamcoeff([float]): Chebyshev polynomial coefficients to the best-fit
            mercury and xe lamp spectrum. lamcoeff should be used as the 
            wavelength solution.
        lamrms(float): The RMS residual for the best wavelength solution.

    Examples:
        You should use the values as follows:

        extractions = np.load(...)

        extraction = extractions[500] # 500 is arbitrary for this example
        pixel = 30
        X = extraction.xrange[0] + pixel
        wavelength = chebval(pixel, extraction.lamcoeff)

        print "The spectrum starting at %i and pixel %i is %f" % \
            (extraction.xrange[0], pixel, wavelength)

    '''

    seg_id = None 
    ok = None
    xrange = None
    yrange = None
    poly = None
    spec = None
    hg_lines = None

    def __init__(self, seg_id=None, ok=None, xrange=None, 
        yrange=None, poly=None, spec=None,
        hg_lines = None):
        

        self.seg_id = seg_id
        self.ok = ok
        self.xrange = xrange
        self.yrange = yrange
        self.poly = poly
        self.spec = spec
        self.hg_lines = hg_lines