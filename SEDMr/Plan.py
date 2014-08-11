
import argparse
import pdb
import numpy as np
import pylab as pl
import pyfits as pf
import sys


import NPK.Fit 
import NPK.Bar as Bar
from astropy.table import Table 

from scipy.spatial import KDTree 
import scipy.signal as SG


from numpy.polynomial.chebyshev import chebfit, chebval

import Extraction
reload(NPK.Fit)
reload(Extraction)


def extract_info(infiles):
    

    headers = []

    update_rate = len(infiles) / Bar.setup()
    for ix, file in enumerate(infiles):
        if ix % update_rate == 0: Bar.update()
        FF = pf.open(file)
        FF[0].header['filename'] = file
        headers.append(FF[0].header)
    
    Bar.done()
    
    return sorted(headers, key=lambda x: x['JD'])

def identify_observations(headers):
    '''Return a list of object name, observation number, and list of files.

    e.g. will return:

    {'STD-BD+25d4655': {1: ['...']}, {2: ['...']}, 
           'PTF14dvo': {1: ['...', '...']}}
    
    where STD-BD+25d4655 was observed at the beginning and end of night. SN
    14dov was observed once with A-B.'''
    JD = 0

    objcnt = {}
    objs = {}

    curr = ""
    for header in headers:
        if header['JD'] < JD:
            raise Exception("Headers not sorted by JD")
        JD = header['JD']

        fname = header['filename']
        obj = header['OBJECT']
        name = header['NAME']
        exptime = header['exptime']
        if obj.startswith('Calib:'): continue
        if "Focus:" in obj: continue
        if "dark" in obj: continue
        if obj.rstrip() == "": continue

        if '[A]' in obj:
            cnt = objcnt.get(name, 0) + 1
            objcnt[name] = cnt
            objs[name] = {cnt: [fname]}
        else:
            cnt = objcnt[name]
            objs[name][cnt].append(fname)


    return objs


make_preamble = '''
EXTSINGLE = @echo
EXTPAIR = @echo
'''

def MF_single(objname, obsnum, file):
    '''Create the MF entry for a observation with a single file. '''

    tp = {'objname': objname, 'obsfile': files[0]}
    if obsnum == 1: tp['num'] = ''
    else: tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp

    return '''%(outname)s: %(obsfile)s
\t$(EXTSINGLE) %(obsfile)s --outfile %(outname)s\n''' %  tp, "%(outname)s" % tp

def MF_AB(objname, obsnum, A, B):
    '''Create the MF entry for an A-B observation'''

    tp = {'objname': objname, 'A': A, 'B': B}
    if obsnum == 1: tp['num'] = ''
    else: tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp
    return '''%(outname)s: %(A)s %(B)s
\t$(EXTPAIR) %(A)s %(B)s --outfile %(outname)s\n''' %  tp, "%(outname)s" % tp



def to_makefile(objs):
    
    MF = ""

    all = "\nall: "
    for objname, observations in objs.iteritems():
        
        objname = objname.replace(" ", "_")
        objname = objname.replace("+", "_")
        objname = objname.replace(")", "_")
        objname = objname.replace("(", "_")
        objname = objname.replace("[", "_")
        objname = objname.replace("]", "_")
        print objname, observations
        for obsnum, obsfiles in observations.iteritems():
            if len(obsfiles) == 1:
                m,a = MF_single(objname, obsnum, obsfiles[0])
                MF += m
                all += a + " "
            if len(obsfiles) == 2:
                m,a = MF_AB(objname, obsnum, obsfiles[0], obsfiles[1])
                MF += m
                all += a + " "
                
    all += "\n\n"


    f = open("Makefile", "w")
    f.write(make_preamble + all + MF)
    f.close()

def make_plan(headers):
    '''Convert headers to a makefile
    
    Assumed headers sorted by JD'''

            
    objs = identify_observations(headers)
    to_makefile(objs)

if __name__ == '__main__':

    files = sys.argv[1:]
    to_process = extract_info(files)

    objs = make_plan(to_process)


