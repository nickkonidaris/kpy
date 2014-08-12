
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
import NPK.Standards as Stds
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
EXTSINGLE =  ~/spy /scr2/npk/PYTHON/SEDMr/Extracter.py 
EXTPAIR =  ~/spy /scr2/npk/PYTHON/SEDMr/Extracter.py 
'''

def MF_single(objname, obsnum, file, standard):
    '''Create the MF entry for a observation with a single file. '''

    print objname, obsnum, file
    tp = {'objname': objname, 'obsfile': file}
    if obsnum == 1: tp['num'] = ''
    else: tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp

    if standard is None: tp['STD'] = ''
    else: tp['STD'] = "--std %s" % (standard)

    first = '''%(outname)s: fine.npy %(obsfile)s 
\t$(EXTSINGLE) fine.npy --A %(obsfile)s --outname %(outname)s %(STD)s --nsighi 1.1\n''' % tp
    second = '''\t$(EXTSINGLE) CORR --A %(outname)s --std %(objname)s --outname CORR.npy\n''' %  tp
    fn = "%(outname)s" % tp

    if standard is None: return first, fn
    else: return first+second, fn

    

def MF_AB(objname, obsnum, A, B):
    '''Create the MF entry for an A-B observation'''

    print objname, obsnum, A, B
    tp = {'objname': objname, 'A': A, 'B': B}
    if obsnum == 1: tp['num'] = ''
    else: tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp
    return '''%(outname)s: fine.npy %(A)s %(B)s
\t$(EXTPAIR) fine.npy %(A)s %(B)s --outname %(outname)s\n''' %  tp, "%(outname)s" % tp



def to_makefile(objs):
    
    MF = ""

    all = "\nall: stds "
    stds = "\nstds: "
    for objname, observations in objs.iteritems():
        
        objname = objname.replace(" ", "_")
        objname = objname.replace(")", "_")
        objname = objname.replace("(", "_")
        objname = objname.replace("[", "_")
        objname = objname.replace("]", "_")
        for obsnum, obsfiles in observations.iteritems():
            if len(obsfiles) == 1:
                standard = None
                if objname.startswith("STD-"):
                    pred = objname[4:].rstrip().lower().replace("+","")
                    print pred
                    if pred in Stds.Standards:
                        standard = pred
                    else: standard = None

                m,a = MF_single(objname, obsnum, obsfiles[0], standard=standard)

                if standard is not None:
                    stds += "%s " % (a)

                MF += m
                all += a + " "
            if len(obsfiles) == 2:
                m,a = MF_AB(objname, obsnum, obsfiles[0], obsfiles[1])
                MF += m
                all += a + " "
                
    all += "\n\n"
    stds += "\n"


    f = open("Makefile", "w")
    f.write(make_preamble + stds + all + MF)
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


