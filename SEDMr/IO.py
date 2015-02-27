
import os
import pyfits as pf
import NPK.Util as UU

def readfits(path):
    ''' Read fits file at path or path.gz '''

    if not os.path.exists(path):
        if os.path.exists("%s.gz" % path):
            path += ".gz"
        else:
            raise Exception("The file at path %s or %s.gz does not exist" % (path, path))

    hdulist = pf.open(path)
    
    return hdulist

def writefits(towrite, fname, no_lossy_compress=False, clobber=False):


    if type(towrite) == pf.PrimaryHDU:
        list = pf.HDUList(towrite)
    elif type(towrite) == pf.HDUList:
        list = towrite
    else:
        list = pf.HDUList(pf.PrimaryHDU(towrite))

    if no_lossy_compress: 
        list.writeto(fname.rstrip(".gz"), clobber=clobber)
        return
    
    list[0].data = UU.floatcompress(list[0].data)
    list.writeto(fname.rstrip(".gz"), clobber=clobber)

    os.system("gzip --force %s" % fname.rstrip(".gz"))
    
