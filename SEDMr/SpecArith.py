
import argparse
import numpy as np
import os
import pylab as pl
import pyfits as pf
from scipy.interpolate import interp1d
import sys


def specDiv(A, B, out):
    ''' Divide spectra in A, B and store to out '''

    s_A = np.load(A)[0]
    s_B = np.load(B)[0]

    t_A, t_B = s_A['exptime'], s_B['exptime']
    if t_A != t_B:
        print "Exposure times do not match (%s v %s). This may be a problem." % (t_A, t_B)

    result = {}
    for k,v in s_A.iteritems():
        if k in ['meta', 'Extinction Correction', 'doc']: continue
        if k not in s_B:
            print "A has key %s but B does not. Skipping." % k
            continue

        if "dlam" in k or "extinction_corr" in k:
            result[k] = (s_A[k] + s_B[k])/2.0
        elif "spectra" in k or "object_spaxel_ids" in k:
            result[k] = np.concatenate((s_A[k], s_B[k]))
        else:
            result[k] = s_A[k] + s_B[k]

    f2 = interp1d(s_B['nm'], s_B['ph_10m_nm'], bounds_error=0, fill_value=np.nan)

    result['ph_10m_nm'] = (s_A['ph_10m_nm'] / f2(s_A['nm']))
    result['nm'] = s_A['nm']

    if s_A.has_key('skyph') and s_B.has_key('skyph') and s_A.has_key('skynm') and s_B.has_key('skynm'):
        s2 = interp1d(s_B['skynm'], s_B['skyph'], bounds_error=0, fill_value=np.nan)
        result['skyph'] = (s_A['skyph'] + s2(s_A['skynm']))/2.0
        result['skynm'] = s_A['skynm']

    if s_A.has_key('var') and s_B.has_key('var'):
        v2 = interp1d(s_B['nm'], s_B['var'], bounds_error=0, fill_value=np.nan)
        result['var'] = (s_A['var'] + v2(s_A['nm']))/2.0

    result['meta_1'] = s_A['meta']
    result['meta_2'] = s_B['meta']
    result['meta'] = "Result of SpecArith."

    np.save(out, [result])

 

def specAdd(A, B, out):
    ''' Add spectra in A, B and store to out '''

    s_A = np.load(A)[0]
    s_B = np.load(B)[0]

    t_A, t_B = s_A['exptime'], s_B['exptime']
    if t_A != t_B:
        print "Exposure times do not match (%s v %s). This may be a problem." % (t_A, t_B)

    result = {}
    for k,v in s_A.iteritems():
        if k in ['meta', 'Extinction Correction', 'doc']: continue
        if k not in s_B:
            raise Exception("A has key %s but B does not. Quitting" % k)

        if "dlam" in k or "extinction_corr" in k:
            result[k] = (s_A[k] + s_B[k])/2.0
        elif "spectra" in k or "object_spaxel_ids" in k:
            result[k] = np.concatenate((s_A[k], s_B[k]))
        else:
            result[k] = s_A[k] + s_B[k]


    print "other"

    f2 = interp1d(s_B['nm'], s_B['ph_10m_nm'], bounds_error=0, fill_value=np.nan)
    s2 = interp1d(s_B['nm'], s_B['skyph'], bounds_error=0, fill_value=np.nan)

    result['ph_10m_nm'] = (s_A['ph_10m_nm'] + f2(s_A['nm']))/2.0
    result['nm'] = s_A['nm']
    result['skyph'] = (s_A['skyph'] + s2(s_A['nm']))/2.0

    if s_A.has_key('var'):
        v2 = interp1d(s_B['nm'], s_B['var'], bounds_error=0, fill_value=np.nan)
        result['var'] = (s_A['var'] + v2(s_A['nm']))/2.0

    result['meta_1'] = s_A['meta']
    result['meta_2'] = s_B['meta']
    result['meta'] = "Result of SpecArith."

    np.save(out, [result])

            

 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''SpecArith.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('--operation', type=str, help='Operation to perform',
        default='+')
    parser.add_argument('A', type=str, help='First term', default="must set A")
    parser.add_argument('B', type=str, help='Second term', default="must set B")
    parser.add_argument('outname', type=str, help='Output name', default="Must set outname")

    args = parser.parse_args()
    err = ""
    if not os.path.isfile(args.A):
        err += "File '%s' does not exist.\n" % args.A
    if not os.path.isfile(args.B):
        err += "File '%s' does not exist.\n" % args.B
    if os.path.isfile(args.outname):
        err += "File '%s' exists. Use a different file name.\n" % args.outname

    if err != "":
        print err
        sys.exit(1)

    print "%s %s %s > %s" % (args.A, args.operation, args.B, args.outname)
    if args.operation == '+':
        specAdd(args.A, args.B, args.outname)
    if args.operation == '/':
        specDiv(args.A, args.B, args.outname)
    else:
        print "%s not recognized as an op." % args.operation
