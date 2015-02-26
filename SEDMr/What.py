
import shutil
import glob
import os
import pyfits as pf

def go(files):

    for file in files:
        FF = pf.open(file)
        header = FF[0].header
        header['FNAME'] = file
        try: print "%(FNAME)28s (%(AIRMASS)1.3f/%(ADCSPEED)1.1f/%(EXPTIME)s s): %(OBJECT)-30s" % (header)
        except:
            print "%28s : ?" % (file)


           
if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        raise Exception("not enough arguments")
        

    go(sys.argv[2:])
        

