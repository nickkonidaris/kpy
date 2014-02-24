import numpy as np
from scipy.optimize import curve_fit

def wavelength_sol(pix, a, b, c, scale, shift):
    '''The wavelength solution function converts pixels to nm'''

    return (a * b**pix + c*pix**2) * scale - shift

def shift_pixels(lam, n_pix):
    """Refits the wavelength fitting function and pixel shifts the solution

    Args:
        [lam]: Array with the wavelength solution in nm.
        n_pix: The amount of pixel shift to apply

    Returns:
        Array [len(lam)] with the new wavelength solution

    """

    xx = np.arange(len(lam))
    pars = [1000., 239./240.0, 0.0, 1.0, 0.0]
    ok = np.isfinite(lam)
    popt, pcov = curve_fit(wavelength_sol,
                            xx[ok],
                            lam[ok],
                            pars)
        
    return wavelength_sol(xx-n_pix, *popt)

