
import os
import sys
import numpy as np
import pyfits as pf
import NPK.Standards as SS
import NPK.Atmosphere as Atm
from pylab import *
from scipy.interpolate import interp1d

stds = ['hz44', 'bd28d' ,'bd33d', 'feige34']
plan = {}
def handle_date(datename):
    ''' Walk directory like/scr2/npk/sedm/reduced/hz44/2014_may_03_20_44_06/v00 for stds'''
    for path, directories, files in os.walk("/scr2/npk/sedm/reduced"):
        if datename not in path: continue

        found_std = False
        for std in stds:
            if std in path: 
                found_std = True
                break
       
        if not found_std: continue
        if len(files) < 3: continue

        # Path will be e.g., /scr2/npk/sedm/reduced/hz44/2014_may_03_20_41_03/v00
        sp = path.split("/")

        name, datetime, ver = sp[-3:]
        try: yy, mnth, dy, hh, mm, ss = datetime.split("_")
        except: continue

        plan[name] = (path, ver)

    return plan
        

def get_std_spec(name):
    '''Returns an interpolating function (wavelength_nm) for std star named name'''

    for stdname, spec in SS.Standards.iteritems():
        if name in stdname:
            return interp1d(spec[:,0]/10, spec[:,1], bounds_error=False)

    return None
            
        

def plot_spectra(plan, date):
    '''take plan[std name] -> (path, ver) and plot'''

    figure(1) ; clf()
    figure(2) ; clf()
    xlim(360,980)

    for stdname, todo in plan.iteritems():

        print stdname,todo
        path,ver = todo
        sp = path.split("/")
        name, datetime, ver = sp[-3:]

        FF = pf.open(os.path.join(path, "%s.fits" % name))
        dat = FF[0].data
        itime = FF[0].header['exptime']
        airmass = FF[0].header['airmass']
        l,s = dat[0,:], dat[2,:]/itime
        magperairmass = Atm.ext(l*10)
        mag_ext = magperairmass * airmass
        ext = 10**(-magperairmass/2.5)
        s /= ext
        sourcefun = interp1d(l, s)
        stdfun = get_std_spec(name)
        figure(1)
        plot(l,s)
        plot(l, stdfun(l))
        figure(2)
        corr = stdfun(l)/s
        semilogy(l, corr)

        try:
            all_corr += stdfun(all_l)/sourcefun(all_l)
        except:
            all_l = l
            all_corr = stdfun(l)/s

    all_corr /= len(plan)
    figure(2)
    legend(plan.keys())

    figure(3)
    ok = np.isfinite(all_corr)
    all_l = all_l[ok]
    all_corr = all_corr[ok]

    ok = (all_l > 370) & (all_l < 920) & (all_corr > 0)
    ff = interp1d(all_l[ok], all_corr[ok],bounds_error = False)
    semilogy(all_l,all_corr,'.')
    semilogy(all_l, ff(all_l))

    np.save('correction_%s' % date, [all_l, all_corr])


    show()


if __name__ == '__main__':
    plan = handle_date(sys.argv[1])
    plot_spectra(plan, sys.argv[1])

