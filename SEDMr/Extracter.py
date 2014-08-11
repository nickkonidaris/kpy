
import sys
import numpy as np
import pylab as pl


from scipy.interpolate import interp1d
import Extraction
reload(Extraction)


def plot_spectra_medians(spectra, low=0, hi=0):


    ms = []

    for spectrum in spectra:
        if spectrum.__dict__.has_key('spec') and spectrum.spec is not None \
            and spectrum.lamcoeff is not None:
            l,s = spectrum.get_flambda()
            ms.append(np.median(s))

    ms = np.array(ms)

    pl.figure(1)
    pl.clf()
    pl.plot(ms,'.')
    pl.axhline(low,color='orange')
    pl.axhline(hi,color='red')
    pl.show()
        

def plot_selected_spectra(spectra, low=0, hi=0):
    
    to_extract =[]
    for spectrum in spectra:
        if spectrum.__dict__.has_key('spec') and spectrum.spec is not None \
            and spectrum.lamcoeff is not None:
            l,s = spectrum.get_flambda()

            m = np.median(s)

            if (m>hi) or (m<low):
                to_extract.append(spectrum)


    print len(to_extract)
    pl.figure(2)
    pl.clf()
    l_grid = None
    s_grid = None
    for spectrum in to_extract:
        l,s = spectrum.get_flambda()
        if np.median(s) < 0: pon = -1.0
        else: pon = 1.0

        if l_grid is None:
            l_grid = l
            s_grid = [s*pon]
        else:
            fun = interp1d(l,s*pon, bounds_error=False,fill_value=0)
            s_grid.append(fun(l_grid))
            
            

    medspec = np.mean(s_grid, 0)
    pl.figure(3)
    pl.clf()
    pl.step(l_grid,medspec)
    yl = pl.ylim()
    pl.show()
    pl.figure(2)
    pl.clf()
    s_grid = np.array(s_grid)
    pl.imshow(s_grid,vmin=yl[0], vmax=yl[1])
    pl.colorbar()
    pl.show()


if __name__ == '__main__':
    
    spectra = sys.argv[1]
    lowval = float(sys.argv[2])
    hival = float(sys.argv[3])
    
    spectra = np.load(spectra)
    plot_spectra_medians(spectra, low=lowval, hi=hival)
    plot_selected_spectra(spectra, low=lowval, hi=hival)


