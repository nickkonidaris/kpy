import argparse
import numpy as np
import pylab as pl
import pyfits as pf
import sys

import astropysics


def imcombine(filelist, out, bpmask=None, reject="none", nlow=None,
        nhigh=None):

    '''Convenience wrapper around IRAF task imcombine

    Args:
        filelist: The list of files to imcombine
        out: The full path to the output file
        bpmask: The full path to the bad pixel mask
        reject: none, minmax, sigclip, avsigclip, pclip
        nlow,nhigh: Parameters for minmax rejection, see iraf docs
    
    Returns:
        None

    Side effects:
        Creates the imcombined file at location `out'
    '''

    #TODO: REMOVE Iraf and use python instead. STSCI Python has
    # A builtin routine.
    from pyraf import iraf
    iraf.images()


    filelist = [("%s[0]" % f) for f in filelist]
    pars = iraf.imcombine.getParList()
    iraf.imcombine.unlearn()

    path = "flatcombine.lst"
    f = open(path, "w")
    for file in filelist:
        f.write(file + "\n")
    f.close()

    s = ("%s," * len(filelist))[0:-1]
    s = s % tuple(filelist)

    if reject == 'minmax':
        t = iraf.imcombine("@%s" % path, out, combine="average",
            reject=reject, nlow=nlow, nhigh=nhigh)
    elif reject == 'sigclip':
        t = iraf.imcombine("@%s" % path, out, combine="average",
            reject=reject, lsigma=nlow, hsigma=nhigh)
    else:
        t = iraf.imcombine(s, out, Stdin=filelist, Stdout=1, combine="average",
            reject=reject)

    iraf.imcombine.setParList(pars)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Imcombine.py performs:

        1) Median combination
        2) Mean combine
        3) Mean combine w/ sigma clipping

    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('--files', type=str, nargs='*', default=[])
    parser.add_argument('--Nhi', type=float, default=None)
    parser.add_argument('--Nlo', type=float, default=None)
    parser.add_argument('--reject', type=str, default='none')
    parser.add_argument('--outname', type=str, default=None)
    args = parser.parse_args()

    filelist = args.files
    out = args.outname
    if args.outname is None:
        print "Set --outname"

    imcombine(filelist, out, bpmask=None, reject=args.reject, 
                nlow=args.Nlo, nhigh=args.Nhi)


