
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
        if 'JD' not in FF[0].header:
            print "Skipping %s" % file
            continue
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
        if "GXN:STOW" in name: continue
        if obj.rstrip() == "": continue
        name= name.replace(" ", "_")
        name= name.replace(")", "_")
        name= name.replace("(", "_")
        name= name.replace("[", "_")
        name= name.replace("]", "_")
        name= name.replace("/", "_")
        name= name.replace(":", "_")

        if '[A]' in obj or name not in objcnt:
            cnt = objcnt.get(name, 0) + 1
            objcnt[name] = cnt
            objs[name] = {cnt: [fname]}
        else:
            try: cnt = objcnt[name]
            except: 
                import pdb
                pdb.set_trace()
            objs[name][cnt].append(fname)


    return objs


make_preamble = '''
EXTSINGLE =  ~/spy /scr2/npk/PYTHON/SEDMr/Extracter.py 
ATM =  ~/spy /scr2/npk/PYTHON/SEDMr/AtmCorr.py 
EXTPAIR =  ~/spy /scr2/npk/PYTHON/SEDMr/Extracter.py 
FLEX = ~/spy /scr2/npk//PYTHON/SEDMr/Flexure.py
'''

def MF_single(objname, obsnum, file, standard=None):
    '''Create the MF entry for a observation with a single file. '''

    #print objname, obsnum, file
    tp = {'objname': objname, 'obsfile': file}
    if obsnum == 1: tp['num'] = ''
    else: tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp

    if standard is None: tp['STD'] = ''
    else: tp['STD'] = "--std %s" % (standard)
    tp['flexname'] = "flex_%s.npy" % (file.rstrip(".fits"))
    first = '''%(outname)s: fine.npy %(obsfile)s %(flexname)s
\t$(EXTSINGLE) fine.npy --A %(obsfile)s --outname %(outname)s %(STD)s --nsighi 1.1 --correction std-correction.npy\n''' % tp
    second = '''corr_%(outname)s: %(outname)s
\t$(ATM) CORR --A %(outname)s --std %(objname)s --outname corr_%(outname)s\n''' %  tp
    fn = "%(outname)s" % tp

    if standard is None: return first, fn
    else: return first+second, fn

    

def MF_AB(objname, obsnum, A, B):
    '''Create the MF entry for an A-B observation'''

    #print objname, obsnum, A, B
    tp = {'objname': objname, 'A': A, 'B': B}
    if obsnum == 1: tp['num'] = ''
    else: tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp
    tp['flexname'] = "flex_%s.npy flex_%s.npy" % (A.rstrip('.fits'), 
        B.rstrip('.fits')) 
    tp['bgdnameA'] = "bgd_%s.npy" % (A.rstrip('.fits'), 
        B.rstrip('.fits')) 

    return '''%(outname)s: fine.npy %(A)s %(B)s %(flexname)s
\t$(EXTPAIR) fine.npy --A %(A)s --B %(B)s --outname %(outname)s --correction std-correction.npy\n''' %  tp, "%(outname)s" % tp


def MF_ABCD(objname, obsnum, files): 
    '''Create the MF entry for an A-B observation'''

    A,B,C,D = files
    tp = {'objname': objname, 'A': A, 'B': B, 'C': C, 'D': D}
    if obsnum == 1: tp['num'] = ''
    else: tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp
    tp['flexname'] = "flex_%s.npy flex_%s.npy flex_%s.npy flex_%s.npy" % (
        A.rstrip('.fits'), 
        B.rstrip('.fits'),
        C.rstrip('.fits'),
        D.rstrip('.fits'))
    return '''%(outname)s: fine.npy %(A)s %(B)s %(C)s %(D)s %(flexname)s
\t$(EXTPAIR) fine.npy --A %(A)s --B %(B)s --C %(C)s --D %(D)s --outname %(outname)s --correction std-correction.npy\n''' %  tp, "%(outname)s" % tp



def to_makefile(objs):
    
    MF = ""

    all = ""
    stds = ""
    
    flexures = ""
    
    for objname, observations in objs.iteritems():
        
        objname = objname.replace(" ", "_")
        objname = objname.replace(")", "_")
        objname = objname.replace("(", "_")
        objname = objname.replace("[", "_")
        objname = objname.replace("]", "_")
        for obsnum, obsfiles in observations.iteritems():
            if objname.startswith("STD-"):
                pred = objname[4:].rstrip().lower().replace("+","")
                print pred, Stds.Standards.keys()
                if pred in Stds.Standards:
                    print "FOUND"
                    standard = pred

                    for ix, obsfile in enumerate(obsfiles):
                        m,a = MF_single(objname, ix+1, obsfile, standard=standard)
                        MF += m
                        all += a + " "
                else: standard = None

            elif len(obsfiles) == 1:
                standard = None

                m,a = MF_single(objname, obsnum, obsfiles[0])

                if standard is not None:
                    stds += "corr_%s " % (a)

                MF += m
                all += a + " "
            elif len(obsfiles) == 2:
                m,a = MF_AB(objname, obsnum, obsfiles[0], obsfiles[1])
                MF += m
                all += a + " "
            elif len(obsfiles) == 4:
                m,a = MF_ABCD(objname, obsnum, obsfiles)
                MF += m
                all += a + " "

            for file in obsfiles:
                on = "flex_" + file.rstrip(".fits") + ".npy"
                flexures += "%s: \n\t$(FLEX) fine.npy %s --outfile %s\n\n"\
                    % (on, file, on)
                
    stds += " "


    f = open("Makefile", "w")
    clean = "\nclean:\n\trm %s %s\n\n" % (all, stds)
    corr = "\nstd-correction.npy: %s\n\t$(ATM) SUM --outname std-correction.npy --files %s\n" % (stds, stds)
    f.write(make_preamble + "stds: %s std-correction.npy\n\n" % (stds) +
        "\nall: stds %s %s" % (all, clean) + "\n\n" +
        corr + MF + "\n\n" + flexures)
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


