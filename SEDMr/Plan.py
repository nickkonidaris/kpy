
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

    update_rate = 1. * len(infiles) / Bar.setup()
    for ix, file in enumerate(infiles):
        if ix % update_rate == 0: Bar.update()
        FF = pf.open(file)
        FF[0].header['filename'] = file
        if 'JD' not in FF[0].header:
            #print "Skipping %s" % file
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
    calibs = {}

    curr = ""
    for header in headers:
        if header['JD'] < JD:
            raise Exception("Headers not sorted by JD")
        JD = header['JD']

        fname = header['filename']
        obj = header['OBJECT']
        name = header['NAME']
        exptime = header['exptime']
        adcspeed = header['ADCSPEED']
        if "test" in obj: continue
        if "Calib" in obj:

            def appendToCalibs(Str):

                if Str in obj:
                    if "bias" in Str:
                        Str = "%s%1.1f" % (Str, adcspeed)
                        prefix = ""
                        suffix = ""
                    else:
                        prefix = "crr_b_"
                        suffix = ""

                    calibs[Str] = calibs.get(Str, [])
                    calibs[Str].append(prefix + fname + suffix)

            appendToCalibs("bias")
            appendToCalibs("dome")
            appendToCalibs("Xe")
            appendToCalibs("Hg")
            appendToCalibs("Cd")
            appendToCalibs("Ne")
            appendToCalibs("twilight")

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
            vals = objs.get(name, {})
            objcnt[name] = cnt
            vals[cnt] = [fname]
            objs[name] = vals
        else:
            try: cnt = objcnt[name]
            except: 
                import pdb
                pdb.set_trace()
            objs[name][cnt].append(fname)


    print "-- Calibrations --"
    for k,v in calibs.iteritems():
        print "%15s : %2.0i" % (k, len(v))

    return objs, calibs


make_preamble = '''


PY = /home/npk/spy
PYC = /scr2/npk/PYTHON/SEDM
EXTSINGLE =  $(PY) $(PYC)r/Extracter.py 
ATM =  ~/spy /scr2/npk/PYTHON/SEDMr/AtmCorr.py 
EXTPAIR =  ~/spy /scr2/npk/PYTHON/SEDMr/Extracter.py 
FLEXCMD = ~/spy /scr2/npk//PYTHON/SEDMr/Flexure.py
DEBIAS = ~/spy /scr2/npk//PYTHON/SEDMr/Debias.py
IMCOMBINE = ~/spy /scr2/npk//PYTHON/SEDMr/Imcombine.py

BSUB = $(PY) $(PYC)/Bias.py
BGDSUB =  $(PY) $(PYC)r/SubtractBackground.py
CRRSUB =  /scr2/npk/Ureka/variants/common/bin/PyCosmic --fwhm=2 --iter 4 --rlim 1.0

SRCS = $(wildcard ifu*fits)
BIAS = $(addprefix b_,$(SRCS))
CRRS = $(addprefix crr_,$(BIAS))
BACK = $(addsuffix .gz,$(addprefix bs_,$(CRRS)))
FLEX = $(subst .fits,.npy,$(addprefix flex_,$(BACK)))
TO_FLATTEN = %(FLATFILES)s
FLATTED = $(addprefix flat_,$(TO_FLATTEN))

bias: bias0.1.fits bias2.0.fits $(BIAS)
bgd: $(BGD) bias
crrs: $(CRRS) 
back: $(BACK)


$(BIAS): bias0.1.fits bias2.0.fits
	$(BSUB) $(subst b_,,$@)

$(CRRS): 
	$(CRRSUB) $(subst crr_,,$@) $@ mask$@ 5.0

$(BACK): 
	$(BGDSUB) fine.npy $(subst .gz,,$(subst bs_,,$@))
    

seg_dome.fits: dome.fits
	~/spy /scr2/npk/PYTHON/SEDMr/SexLamps.py dome.fits

seg_Hg.fits: Hg.fits
	~/spy /scr2/npk/PYTHON/SEDMr/SexSpectra.py Hg.fits

dome.fits_segments.npy: seg_dome.fits
	~/spy /scr2/npk/PYTHON/SEDMr/FindSpectra.py seg_dome.fits dome.fits dome.fits_segments

rough.npy: dome.fits_segments.npy seg_Hg.fits
	~/spy /scr2/npk/PYTHON/SEDMr/Wavelength.py rough --hgfits Hg.fits --hgcat cat_Hg.fits.txt --dome dome.fits_segments.npy --outname rough 

fine.npy: rough.npy
	~/spy /scr2/npk/PYTHON/SEDMr/Wavelength.py fine --xefits Xe.fits --hgfits Hg.fits --hgassoc assoc_Hg.npy --outname fine

cube.npy: fine.npy
	~/spy /scr2/npk/PYTHON/SEDMr/Cube.py fine.npy --step make --outname cube.npy

bs_twilight.fits.gz: twilight.fits fine.npy
	$(BGDSUB) fine.npy twilight.fits

bs_dome.fits.gz: dome.fits fine.npy
	$(BGDSUB) fine.npy dome.fits

flat.fits: bs_twilight.fits.gz bs_dome.fits.gz
	~/spy /scr2/npk/PYTHON/SEDMr/MakeFlat.py bs_twilight.fits.gz --toflat2 bs_dome.fits.gz

wave: fine.npy
cube: cube.npy

# TODO: Change the dependencies to be more sophisticated.
# It is possible to add rules to both flex and flat so that only
# the necessary files are created
flex: back $(FLEX)
flat: back flat.fits $(FLATTED)

$(FLEX): cube.npy
	$(eval OUTNAME = $(subst .gz,,$@))
	$(FLEXCMD) cube.npy $(subst flex_,,$(subst npy,fits,$@)) --outfile $(OUTNAME)

$(FLATTED): flat.fits
	$(eval LSRC = $(subst flat_,,$@))
	$(eval LFLEX = $(subst fits.gz,npy,$(addprefix flex_,$(LSRC))))
	~/spy /scr2/npk/PYTHON/SEDMr/DivideFlat.py flat.fits $(LSRC) --flexnpy $(LFLEX) --outfile $@

'''

def MF_imcombine(objname, files, dependencies=""):
    
    filelist = " ".join(["%s " % file for file in files])
    first = "%s.fits: %s %s\n" % (objname, filelist, dependencies)

    if len(filelist) > 7:
        reject = "sigclip"
    else:
        reject = "none"
    second = "\t$(IMCOMBINE) --outname %s.fits --reject %s --Nlo 3 --Nhi 3 --files %s\n" % (objname, reject, filelist)

    return  first+second+"\n"


def MF_single(objname, obsnum, file, standard=None):
    '''Create the MF entry for a observation with a single file. '''

    tp = {'objname': objname, 'obsfile': "bs_crr_b_%s" % file}
    tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp

    if standard is None: tp['STD'] = ''
    else: tp['STD'] = "--std %s" % (standard)
    tp['flexname'] = "flex_bs_crr_b_%s.npy" % (file.rstrip(".fits"))
    first = '''# %(outname)s\n%(outname)s: cube.npy %(flexname)s %(obsfile)s.gz
\t$(EXTSINGLE) cube.npy --A %(obsfile)s.gz --outname %(outname)s %(STD)s --nsighi 1.1 --correction std-correction.npy

cube_%(outname)s.fits: %(outname)s
\t~/spy /scr2/npk/PYTHON/SEDMr/Cube.py %(outname)s --step extract --outname cube_%(outname)s.fits
''' % tp
    second = '''corr_%(outname)s: %(outname)s
\t$(ATM) CORR --A %(outname)s --std %(objname)s --outname corr_%(outname)s\n''' %  tp
    fn = "%(outname)s" % tp

    if standard is None: return first+"\n", fn
    else: return first+second+"\n", fn 

    

def MF_AB(objname, obsnum, A, B):
    '''Create the MF entry for an A-B observation'''

    #print objname, obsnum, A, B
    tp = {'objname': objname, 'A': A, 'B': B}
    if obsnum == 1: tp['num'] = ''
    else: tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp
    tp['flexname'] = "flex_%s.npy flex_%s.npy" % (A.rstrip('.fits'), 
        B.rstrip('.fits')) 
    tp['bgdnameA'] = "bgd_%s.npy" % (A.rstrip('.fits'))
    tp['bgdnameB'] = "bgd_%s.npy" % (B.rstrip('.fits'))  

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



def to_makefile(objs, calibs):
    
    MF = ""

    all = ""
    stds = ""
    
    flexures = ""

    for calibname, files in calibs.iteritems():
        
        if "bias" not in calibname:
            pass
        MF += MF_imcombine(calibname, files)
        all += "%s.fits " % calibname
    
    flatfiles = []
    for objname, observations in objs.iteritems():

        objname = objname.replace(" ", "_")
        objname = objname.replace(")", "_")
        objname = objname.replace("(", "_")
        objname = objname.replace("[", "_")
        objname = objname.replace("]", "_")

        for obsnum, obsfiles in observations.iteritems():
            flatfiles.append(obsfiles)

            if objname.startswith("STD-"):
                pred = objname[4:].rstrip().lower().replace("+","").replace("-","_")
                if pred in Stds.Standards:
                    standard = pred

                    for ix, obsfile in enumerate(obsfiles):
                        m,a = MF_single(objname, obsnum, obsfile, 
                            standard=standard)
                        MF += m
                        all += a + " "

                else: standard = None
                continue

            for obsnum, obsfile in enumerate(obsfiles):
                standard = None

                m,a = MF_single(objname, obsnum, obsfile)

                if standard is not None:
                    stds += "corr_%s " % (a)

                MF += m
                all += a + " "
            '''
            elif len(obsfiles) == 2:
                m,a = MF_AB(objname, obsnum, obsfiles[0], obsfiles[1])
                MF += m
                all += a + " "
            elif len(obsfiles) == 4:
                m,a = MF_ABCD(objname, obsnum, obsfiles)
                MF += m
                all += a + " "
            '''

    stds += " "


    flattened_flatfiles = sum(flatfiles, [])
    flattened_flatfiles = ["bs_crr_b_%s.gz" % fn for fn in flattened_flatfiles]
    flatfilestr = ' '.join([fn for fn in flattened_flatfiles])

    preamble = make_preamble % {'FLATFILES': flatfilestr}

    f = open("Makefile", "w")
    clean = "\nclean:\n\trm %s %s\n\n" % (all, stds)
    corr = "\nstd-correction.npy: %s\n\t$(ATM) SUM --outname std-correction.npy --files %s\n" % (stds, stds)
    f.write(preamble + "stds: %s std-correction.npy\n\n" % (stds) +
        "\nall: stds %s %s" % (all, clean) + "\n\n" +
        corr + MF + "\n\n" + flexures)
    f.close()

def make_plan(headers):
    '''Convert headers to a makefile
    
    Assumed headers sorted by JD'''

            
    objs, calibs = identify_observations(headers)
    to_makefile(objs, calibs)

if __name__ == '__main__':

    files = sys.argv[1:]
    to_process = extract_info(files)

    objs = make_plan(to_process)


