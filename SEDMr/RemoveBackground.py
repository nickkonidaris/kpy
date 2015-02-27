import argparse
import numpy as np
import tempfile
import os
import pyfits as pf


# catalog_name, output_name
sex_params = \
'''
CATALOG_NAME     {catalog_name} # name of the output catalog
CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
PARAMETERS_NAME  /tmp/sex.sex.param  # name of the file containing catalog contents
 
#------------------------------- Extraction ----------------------------------
 
DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)
DETECT_MINAREA   9              # minimum number of pixels above threshold
DETECT_THRESH    1.0            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH  1.0            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
 
FILTER           N              # apply filter for detection (Y or N)?
FILTER_NAME      default.conv   # name of the file containing the filter
 
                                # NPK: These parameters are tuned such that each
                                # peak is nearly object
DEBLEND_NTHRESH  32             # Number of deblending sub-thresholds
DEBLEND_MINCONT  0.1          # Minimum contrast parameter for deblending
 
CLEAN            Y              # Clean spurious detections? (Y or N)?
CLEAN_PARAM      1.0            # Cleaning efficiency
 
MASK_TYPE        CORRECT        # type of detection MASKing: can be one of
                                # NONE, BLANK or CORRECT
 
#------------------------------ Photometry -----------------------------------
 
PHOT_APERTURES   32              # MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,
                                # <min_radius>
 
SATUR_LEVEL      50000.0        # level (in ADUs) at which arises saturation
SATUR_KEY        SATURATE       # keyword for saturation level (in ADUs)
 
MAG_ZEROPOINT    0.0            # magnitude zero-point
MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)
GAIN             0.0            # detector gain in e-/ADU
GAIN_KEY         GAIN           # keyword for detector gain in e-/ADU
PIXEL_SCALE      1.0            # size of pixel in arcsec (0=use FITS WCS info)
                                # PIXEL_SCALE set by npk
 
#------------------------- Star/Galaxy Separation ----------------------------
 
SEEING_FWHM      2.0            # stellar FWHM in arcsec. Set to 2.5 by NPK 
STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename
 
#------------------------------ Background -----------------------------------
 
BACK_SIZE        200,100 # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE  1              # Background filter: <size> or <width>,<height>
 
BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL
 
#------------------------------ Check Image ----------------------------------
 
CHECKIMAGE_TYPE  -BACKGROUND  SEGMENTATION BACKGROUND  FILTERED
                                # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                # or APERTURES
CHECKIMAGE_NAME  s_{output_name} seg_{output_name} back_{output_name} filtered_{output_name}
 
#--------------------- Memory (change with caution!) -------------------------
 
MEMORY_OBJSTACK  3000           # number of objects in stack
MEMORY_PIXSTACK  300000         # number of pixels in stack
MEMORY_BUFSIZE   1024           # number of lines in buffer
 
#----------------------------- Miscellaneous ---------------------------------
 
VERBOSE_TYPE     QUIET # can be QUIET, NORMAL or FULL
WRITE_XML        N              # Write XML file (Y/N)?
XML_NAME         sex.xml        # Filename for XML output

'''


def handle_path(path):
    name= os.path.basename(path)

    c = sex_params.format(**{"catalog_name": "cat_%s.txt" % name, 
        "output_name": name})

    conf_file = open("/tmp/sedm_sex_conf.sex","w")
    conf_file.write(c)
    conf_file.close()
    os.system("sex -c /tmp/sedm_sex_conf.sex {0}".format(path))


def remove(fname):
    try: os.remove(fname)
    except: pass

def grow(img):

    img[img > 0] = 1
    cp = img[:,:]

    for ix in xrange(img.shape[0]):
        for jx in xrange(img.shape[1]):
            if img[ix,jx] == 1:
                try:
                    cp[ix,jx+1] = 1
                    cp[ix,jx-1] = 1
                    cp[ix+1,jx+1] = 1
                    cp[ix+1,jx-1] = 1
                    cp[ix-1,jx+1] = 1
                    cp[ix-1,jx-1] = 1
                except IndexError:
                    continue

    return cp

def go(paths):

    f = open('/tmp/sex.sex.param', 'w')
    f.write("NUMBER\nX_IMAGE\nY_IMAGE\nFLUX_ISO\nISOAREA_IMAGE\n")
    f.close()

    f = open("default.conv", "w")
    f.write("""CONV NORM
# 2 pix fwhm
1 2 1
2 4 2
1 2 1
""")
    f.close()

    for path in paths:
        handle_path(path)
        # Read background image
        name= os.path.basename(path)
        Bi = pf.open("back_%s" % name)
        Si = pf.open("seg_%s" % name)
        Di = pf.open("%s" % name)

        mask = grow(Si[0].data)
        mask = grow(mask)
        tofill = mask > 0
        if tofill.any():
            Di[0].data[tofill] = Bi[0].data[tofill]
            #Di[0].data[tofill] = np.nan
            Di.writeto("backfill_%s" % name, clobber=True)

        path = path.replace(name, "backfill_%s" % name)
        handle_path(path)

        Di = pf.open("%s" % name)
        Bi = pf.open("back_backfill_%s" % name)

        Di[0].data -= Bi[0].data
        Di.writeto("bs_%s" % name, clobber=True)

        remove("backfill_%s" % name)
        remove("s_backfill_%s" % name)
        remove("filtered_backfill_%s" % name)
        remove("seg_backfill_%s" % name)
        remove("cat_backfill_%s.txt" % name)

        remove("back_%s" % name)
        remove("seg_%s" % name)
        remove("filtered_%s" % name)
        remove("seg__%s" % name)
        remove("cat_%s.txt" % name)


            

        



if __name__ == '__main__':
    
    import sys

    go(sys.argv[1:])
    

