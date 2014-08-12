'''
SegmentationMap class wraps around several functions in SEDM package and
provides a convenient way to handle SEDM data.
'''

import NPK.Atmosphere as Atm
import NPK.PlotHelp as PH
import datetime
import numpy as np
import json
import os
import pyfits
import scipy.io
import matplotlib.pyplot as pl
from matplotlib.backend_bases import KeyEvent 
from matplotlib.backend_bases import PickEvent
from scipy.interpolate import interp1d

import Disp
import SegMap
reload(SegMap)
from SegMap import SegmentationMap


import Extract
from NPK.Standards import Standards

reload(Extract)
reload(Disp)

def clean_header(header):

    todel = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
        'NAXIS0', 'EXTEND', 'BZERO', 'BSCALE']
    for d in todel: 
        try: del header[d]
        except: pass

    return header

def merge_headers(headers):
    h = headers[0]

    i = 1
    for new in headers[1:]:
        for key in new.keys():
            if 'NAXIS' in key: continue
            if key in h:
                if h[key] != new[key]:
                    try: h[key + "_%2.2i" % i] = new[key]
                    except: pass
            else:
                try: h[key + "_%2.2i" % i] = new[key]
                except: pass

    return h

class PositionPicker(object):
    '''Shows the IFU field and allows the user to select an object.

    Takes:
        A plan list that is composed of dictionary that includes:
            infiles [N_spec]: String of paths to .mat files
            name: String of object name
        Optional:
            outdir: path to write the results to, otherwise defaults to
                oudri/[name]/[obsdate]/version#/....
            shifts [N_spec]: Float of pixels to shift spectrum for flexure

    '''

    index = 0   # index into the file list
    spx_ix = 0 # index into spectrum list
    plan = None
    subtract_sky = False
    picked = None # Spectrum position pixed in IFU coordinates
    olines  = [372.7, 486.1, 500.7, 656.3]
    tlines = [761.5, 589.0, 557.7, 435.8, 519.9, 630.0]
    positions=None
    qecurve = None
    output_head = "/scr2/npk/sedm/reduced/"
    show_calib = True
    norm = None
    

    def __init__(self,plan,positions=('OnSkyX','OnSkyY'),
        qefunction=None, stdpaths=None, normpath=None):
        print "Starting picker GUI"
        self.plan = plan
        self.positions=positions
        self.stdpaths = stdpaths
        self.load_stds()
        if normpath is not None:
            self.norm = np.load(normpath)
            self.norm[self.norm < .5] = np.nan
            self.norm[self.norm > 2] = np.nan

        self.fig = pl.figure(1)
        self.fig.canvas.mpl_connect("key_press_event", self)
        self.fig.canvas.mpl_connect("pick_event", self)
        self.fig.canvas.mpl_connect("button_press_event", self)

        self.fig2 = pl.figure(2,figsize=(16,4.5))
        self.fig2.canvas.mpl_connect("key_press_event", self)
        self.fig2.canvas.mpl_connect("pick_event", self)

        self.fig3 = pl.figure(3, figsize=(12,4))
        self.fig3.canvas.mpl_connect("key_press_event", self)
        self.fig3.canvas.mpl_connect("pick_event", self)

        self.fig4 = pl.figure(4, figsize=(8,2.25))
        self.fig4.canvas.mpl_connect("key_press_event", self)
        self.fig4.canvas.mpl_connect("pick_event", self)

        if qefunction is not None:
            self.qecurve = qefunction
        self.index = 0

        self.check_plan()
        self.load()
        self.draw()

    def check_plan(self):
        '''Checks the plan to see if it makes sense and the files load'''

        for el in self.plan:
            print el['name']
            for fname in el['infiles']:
                if not os.path.exists(fname):
                    raise Exception("%s: file does not exist" % fname)

            if 'object_diam' in el:
                if (0 > el['object_diam'] < 300):
                    raise Exception("%s: not appropriate object diam" %
                        el['object_diam'])

            if 'sky_annulus' in el:
                if len(el['sky_annulus']) != 2:
                    raise Exception("%s: not appropriate sky annulus" %
                        el["sky_annulus"])
                


    def dump(self):
        """Write status to file"""
        plan = self.plan[self.index]

        dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if self.status == 'standard':
            outfile = os.path.join(self.outdir,
                "STD-%s.json" % plan['name'])
        else:
            outfile = os.path.join(self.outdir,
                "%s.json" % plan['name'])

        print "Writing to: %s" % outfile
        str = json.dumps({"outfile": outfile,
                        "infiles": plan['infiles'],
                        "object_diam": plan['object_diam'],
                        "sky_annulus": plan['sky_annulus'],
                        "status": self.status,
                        "picked": self.picked,
                        "pixel_shift": self.pixel_shift,
                        "when": dt}, indent=4)

        try:
            f=open(outfile, "w")
            f.write(str)
            f.close()
        except Exception as e:
            raise("Could not write %s: %s" % (outfile, e))

        outfile = os.path.join(self.outdir,
            "%s_all" % plan['name'])
        np.savez(outfile, self.all_obj_spec)

        outfile = os.path.join(self.outdir,
            "%s_sky" % plan['name'])
        np.savez(outfile, self.all_sky_spec)

        
        result = np.array([self.sky_spec[0], self.sky_spec[1],
                            self.obj_spec[1]-self.sky_spec[1]])
        self.RESULTS[self.spx_ix] = result

        header = merge_headers(self.headers)
        pf = pyfits.PrimaryHDU(result,header=header)
        name = self.plan[self.index]['name']
        pf.header["OBJECT"] = name
        pf.header["SPEC"] = self.plan[self.index]['infiles'][self.spx_ix]
        outpath = os.path.join(self.outdir, "%s_%i.fits" % (name,
                self.spx_ix))

        try: os.remove(outpath)
        except: pass
            
        pf.writeto(outpath)

        self.draw_res()

    def draw_res(self):
        ''' Draw the resulting spectrum'''
        pl.figure(3)
        pl.clf()
        pl.xlim([350, 950])
        pl.ylim(-1000,4000)

        allwav = allsky = allobj = None
        headers = []
        for i in xrange(len(self.RESULTS)):
            res = self.RESULTS[i]
            if res is None: continue

            header = self.headers[i]
            header = clean_header(header)
            headers.append(header)
            airmass = header['airmass']

            wave, sky, obj = res.copy()
            skyf = interp1d(wave,sky,fill_value=np.nan,bounds_error=False)
            objf = interp1d(wave,obj,fill_value=np.nan,bounds_error=False)


            if self.qecurve is None: correction = 1.0
            else: 
                print "Applying correction"
                ext = 10**(-Atm.ext(wave*10)*airmass/2.5)
                correction = 1/self.qecurve(wave*10)*ext
                correction=1.

            pl.step(wave,obj*correction)

            if allwav is None: 
                allwav=wave[:]
                allsky=sky[:]
                allobj=obj[:]
            else:
                allsky += skyf(allwav)
                allobj += objf(allwav)

        if self.qecurve is None: correction = 1.0
        else: correction = 1/self.qecurve(allwav*10)

        pl.step(allwav, allobj*correction, linewidth=3)
        #pl.step(allwav, allobj, linewidth=3)
        

        # Raw data
        result = np.array([allwav, allsky, allobj])
        pf = pyfits.PrimaryHDU(result, header=header)
        name = self.plan[self.index]['name']
        pf.header["REDNAME"] = name
        outpath = os.path.join(self.outdir, "%s.fits" % (name))

        try: os.remove(outpath)
        except: pass

        pf.writeto(outpath)
        # Corrected
        if self.qecurve is not None:
            result = np.array([allwav, allsky*correction, allobj*correction])
            pf = pyfits.PrimaryHDU(result, header=header)
            name = self.plan[self.index]['name']
            pf.header["REDNAME"] = name
            outpath = os.path.join(self.outdir, "%s_corr.fits" % (name))
            try: os.remove(outpath)
            except: pass
            pf.writeto(outpath)


    def load(self):
        """Load the Segmentation map"""

        self.spx_ix = 0

        pl.figure(1)
        pl.clf()
        if self.index < 0: self.index =0
        if self.index >= len(self.plan): 
            self.index = len(self.plan)-1
            print "REACHED THE END"

        cur_plan = self.plan[self.index]

        nfiles = len(cur_plan['infiles'])

        self.SM = []
        self.headers = []
        for i in xrange(nfiles):
            self.SM.append(SegmentationMap(cur_plan['infiles'][i],
                positions=self.positions, norm=self.norm))
            
            fits = cur_plan['infiles'][i].replace("shrunk_","").rstrip('_SI.mat')
            FF = pyfits.open(fits)
            self.headers.append(FF[0].header)

        self.status = ["unknown"] * len(self.SM)
        self.picked = [None] * len(self.SM)
        self.pixel_shift = [0.0] * len(self.SM)
        self.RESULTS = [None] * len(self.SM)

        
        fn = cur_plan['infiles'][0]
        fname = fn.split("/")[-1]
        print fname
        if "shrunk" in fname: fname = fname.replace("shrunk","")
        if "_crr" in fname: fname = fname.replace("_crr","")
        if "_s" in fname: fname = fname.replace("_s","")
        if "_b" in fname: fname = fname.replace("_b","")
        if "b_" in fname: fname = fname.replace("b_","")
        if "ifu" in fname: fname = fname.replace("_ifu","")
        fname = fname.rstrip(".fits_SI.mat")

        print fname
        fn = fname

        y,mon,d = (fn[0:4], fn[4:6], fn[6:8])
        h,min,s = (fn[9:11], fn[12:14], fn[15:17])

        month = ["none","jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
            "sep", "oct", "nov", "dec"][int(mon)]

        outprefix = os.path.join(self.output_head, 
            cur_plan['name'], 
            '%s_%s_%s_%s_%s_%s' % (y, month, d, h, min, s))

        try:
            os.makedirs(outprefix)
        except:
            pass

        
        for i in xrange(99):
            outdir = os.path.join(outprefix, "v%2.2i" % i)
            if not os.path.exists(outdir):
                break
            if len(os.listdir(outdir)) == 0:
                break

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        self.outdir = outdir
                
        print ("Loaded. Outdir: %s" % (self.outdir))

    def create_spectrum(self):
        
        plan = self.plan[self.index]
        if self.picked[self.spx_ix] is None: return
        self.SM[self.spx_ix].pixel_shift = self.pixel_shift[self.spx_ix]

        X,Y = self.picked[self.spx_ix]
        header = self.headers[self.spx_ix]
        header = clean_header(header)

        exptime = header['exptime']
        if "object_diam" in plan: D = plan['object_diam']
        else: 
            D = 2
            plan['object_diam'] = D
        header['extdiam'] = (D, 'Extraction diameter')

        if "sky_annulus" in plan: s1,s2 = plan['sky_annulus']
        else: 
            s1,s2 = 4,6
            plan['sky_annulus'] = [s1,s2]
        
        header['sky1'] = (s1, 'Sky inner annulus diameter')
        header['sky2'] = (s2, 'Sky outer annulus diameter')
        header['ext_pos'] = ("%s,%s" % (X,Y), 'Extraction position')
        print "Extracting at: %s,%s: diam %s, sky %s/%s" % (X,Y,D,s1,s2)

        sky_spec, all_sky_spec = self.SM[self.spx_ix].spectrum_in_annulus(
            X,Y,small=s1,
            large=s2)
        wave = sky_spec["wave_nm"]
        obj_spec, all_obj_spec = self.SM[self.spx_ix].spectrum_near_position(X,Y, 
            distance=D, onto=wave)

        sky = [wave, sky_spec["spec_adu"]/len(all_sky_spec)/exptime*len(all_obj_spec)]
        object = [wave, obj_spec["spec_adu"]/exptime]
        self.sky_spec = sky
        self.obj_spec = object
        self.all_sky_spec = all_sky_spec
        self.all_obj_spec = all_obj_spec


    def draw_spectrum(self):
        ''' Draw the spectrum in figure(2) 
        '''

        pl.figure(2)
        pl.clf()
        pl.xlim(350,920)
        if self.picked[self.spx_ix] is None: 
            return

        pl.xlabel("wavelength [nm]")

        obj_spec = self.obj_spec[:]
        sky_spec = self.sky_spec[:]
        wave = obj_spec[0]

        pl.step(obj_spec[0], obj_spec[1], linewidth=1)
        pl.step(sky_spec[0], sky_spec[1], linewidth=1)
        pl.step(sky_spec[0], obj_spec[1]-sky_spec[1], linewidth=2)
        obj_spec[1] -= sky_spec[1]

        correction = 1.0
        if self.qecurve is not None:
            correction = 1/self.qecurve(wave*10)
            correction /= np.median(correction)
            bad = (wave < 400) | (wave > 920) | (correction < 0)
            correction[bad] =np.nan

        if self.show_calib:
            pl.step(wave, obj_spec[1]*correction)
            pl.step(wave, sky_spec[1]*correction)
                
        PH.transparent_legend(['o', 's', 'o-s', 'c x (o-s)', 'c x s'])
        
        pl.xlim(350,920)

        for line in self.olines:
            pl.axvline(line)
        for line in self.tlines:
            pl.axvline(line, color='r')

    def handle_shift(self, xdata, ydata):

        if (xdata < 360) or (xdata > 1000): return

        lines = np.concatenate((self.olines, self.tlines))

        delts = (lines - xdata)
        ix = np.argmin(np.abs(delts))
        print "Closest to %f" % lines[ix]
        
        line = lines[ix]
        delt = delts[ix]
        wave = self.sky_spec[0]

        wix = np.nanargmin(np.abs(wave-line))
        dw = wave[wix]-wave[wix-1]
        print "Delt: {0}, dw: {1}".format(delt, dw)
        self.pixel_shift[self.spx_ix] += delt/dw

        print "pixel shift is: {0}".format(self.pixel_shift[self.spx_ix])
        self.draw_spectrum()

    def __call__(self, event):
        '''Event call handler for Picker gui.'''
        
        if event.name == 'button_press_event':
            print event.xdata, event.ydata
            self.picked[self.spx_ix] = (event.xdata, event.ydata)
            self.create_spectrum()
            self.draw_spectrum()
            self.draw_selection_circle()
            
        elif event.name == 'key_press_event':
            if event.key == '\\':
                print "Shifting"
                self.handle_shift(event.xdata, event.ydata)

            if event.key == '.':
                self.next_spec()

            if event.key == ',':
                self.spx_ix -= 1
                if self.spx_ix < 0: 
                    self.spx_ix = 0
                else:
                    self.draw()
                    self.draw_spectrum()

            if event.key == 'c':
                print "toggle show calib"
                self.show_calib = not self.show_calib
                self.draw_spectrum()

            if event.key == 'n': 
                print "next"
                self.index += 1
                self.load()
                self.draw()
            if event.key == 'p': 
                print "prev"
                self.index -= 1
                self.load()
                self.draw()
            if event.key == '-':
                self.subtract_sky = not self.subtract_sky
                print "Substract sky: %s" % self.subtract_sky
                self.draw()
            if event.key == "u":
                self.status[self.spx_ix] = "unsure"
                self.dump()
                self.draw()
            if event.key == "b":
                self.status[self.spx_ix] = "bad"
                self.dump()
                self.draw()
            if event.key == "o":
                self.status[self.spx_ix] = "ok"
                self.dump()
                self.next_spec()
            if event.key == "s":
                self.status[self.spx_ix] = "standard"
                self.dump()
                self.next_spec()

            if event.key == 'h':
                print """Help---
n - next
p - prev
- - subtract sky
u - unsure: there are targets visible, not sure which is correct.
b - bad: nothing visible
o - ok: target visible
"""
            print event.key
            
        

    def draw(self):
        if self.subtract_sky:
            sky = self.SM[self.spx_ix].sky_median()
            sky_spec = sky['wave_nm'], sky['spec_adu']
        else:
            sky_spec = None

        x,y,v = Extract.segmap_to_img(self.SM[self.spx_ix].SegMap, 
            sky_spec=sky_spec,
            minl=500, maxl=700,positions=self.positions)
        self.Xs = x
        self.Ys = y
        self.Values = v
        self.draw_selection_circle()

    def load_stds(self):
        
        if self.stdpaths is None:
            return

        STDS = []

        for stdpath in self.stdpaths:
            try:
                print stdpath
                FF = pyfits.open(stdpath)
            except:
                print "Ignoring %s" % stdpath
                continue

            name = FF[0].header['OBJECT']
            pl.figure(4)
            pl.clf()
            for std in Standards.keys():
                if name.lower() in std.lower():
                    print name, std
                    ang = Standards[std][:,0]
                    stdflux = Standards[std][:,1] * 1e-16 # to erg/s/cm2/A

                    stdf = interp1d(ang, stdflux, bounds_error = False)

                    # stdf in ADU
                    obsl, obss, obsf = FF[0].data[0,:], FF[0].data[1,:], FF[0].data[2,:]
                    exptime = FF[0].header['EXPTIME']
                    std = (obsf - obss)/exptime
                    print "Std exptime: %s" % exptime

                    to_fit = std/stdf(obsl*10)

                    zero = (obsl> 990)
                    to_fit /= np.median(to_fit)
                    to_fit[zero] = 0.0
                    ok = (np.isfinite(to_fit)) & (obsl<1100) & (obsl>350)
                    poly= np.polyfit(obsl[ok]*10, to_fit[ok], 25)
                    qef = np.poly1d(poly)
                    print poly


                    pl.xlim([350, 950])
                    pl.plot(obsl[ok], to_fit[ok],'.')
                    pl.plot(obsl[ok], qef(obsl[ok]*10))

            pl.xlim([350, 950])


        def qe_clean_functor(qef):
            def to_return(lam):
                v = qef(lam)
                v[lam<3750] = np.nan
                v[lam>9500] = np.nan
                return v
            return to_return

        self.qecurve = qe_clean_functor(qef)

    def draw_selection_circle(self):
        pl.figure(1, figsize=(9,8))

        pl.ion()
        pl.clf()
    
        name = self.plan[self.index]['name']

        if self.positions[0] == 'OnSkyX': diam = 0.5
        else: diam = 1


        pl.title("{0}:{1}. {2} of {3}".format(name,
            self.status[self.spx_ix],
            self.spx_ix,
            len(self.SM)-1))
        if self.subtract_sky:
            vals = self.Values[:]
            cut=25
            vals[vals<-cut/4] = -cut/4
            vals[vals>cut] = cut
        else: vals = self.Values

        v = vals[:]
        #v[v<-100] = -100
        #v[v>50000] = 50000
        pl.scatter(self.Xs,
                    self.Ys,
                    c=vals,
                    s=40,
                    picker=diam,
                    marker='h')
        pl.xlim(-100,2048+100)
        pl.ylim(-100,2048+100)
        pl.colorbar()

        if self.picked[self.spx_ix] is not None:
            X,Y = self.picked[self.spx_ix]
            s1,s2 = self.plan[self.index]['sky_annulus']
            print("Adding circle at {0}/{1}".format(X,Y))
            obj_radius = self.plan[self.index]['object_diam']/2.0
            obj_circle = pl.Circle((X,Y), radius=obj_radius, fill=False,
                color='black',linewidth=2)
            s1_circle = pl.Circle((X,Y), radius=s1/2.0, fill=False,
                color='black',linewidth=2)
            s2_circle = pl.Circle((X,Y), radius=s2/2.0, fill=False,
                color='black',linewidth=2)


            fig = pl.figure(1)
            ax = fig.add_subplot(1,1,1)
            ax.add_patch(obj_circle)
            ax.add_patch(s1_circle)
            ax.add_patch(s2_circle)


    def next_spec(self):
        self.spx_ix += 1
        max_spx_ix = len(self.SM)
        if self.spx_ix >= max_spx_ix: 
            self.spx_ix = max_spx_ix-1
        else:
            self.draw()
            self.draw_spectrum()
            self.draw_selection_circle()



