import numpy as np
import stsci.tools.nmpfit as mpfit

from numpy.polynomial.chebyshev import Chebyshev as CC

def gaussian4(p, x):
    ''' gaussian model
    p[0] -- scale factor
    p[1] -- centroid
    p[2] -- sigma
    p[3] -- offset

    Area: (p[0]-p[3]) x p[2] x np.sqrt(2 x pi)
    '''

    u = (x - p[1])/p[2]
    return p[0]*np.exp(-0.5*u*u) + p[3] 

def gaussian5(p, x):
    ''' gaussian model
    p[0] -- scale factor
    p[1] -- centroid
    p[2] -- sigma
    p[3] -- offset
    p[4] -- slope

    Area: (p[0]-p[3]) x p[2] x np.sqrt(2 x pi)
    '''

    u = (x - p[1])/p[2]
    return p[0]*np.exp(-0.5*u*u) + p[3] + p[4]*x

def sedm_wavelen(p, x):
    ''' SED Machine wavelength function
    '''
    

    A,B,C,D = p

    return 
    return A + B*x + C*x**2 + D*x**3
    return A*B**x  + C*(x-D)**1 


def mpfit_residuals(modelfun):
    '''Returns a residual function for mpfit code'''

    def fun(param, fjac=None, x=None, y=None, error=None):
        '''Generic function'''
        model = modelfun(param, x)
        status = 0

        if error is None:
            return [status, y-model]

        return [status, (y-model)/error]

    return fun


def mpfit_do(residual_fun, # function returned from mpfit_residuals() above
        x, # input x
        y, # input y = f(x)
        parinfo, # initial parameter guess
        error=None,
        quiet=1,
        maxiter=20):
    '''Returns mpfit fit sturcture for residual fun
    
    Args:
        residual_fun: residual_fun from mpfit_residuals
        x, y: x and y data values
        parinfo: Structure containing mpfit pars see
            help(mpfit)

    Example:
        g4res = Fit.mpfit_residuals(Fit.gaussian4)
        parguess = [{'value': 1600}, {'value': 0}, {'value': 2}, {'value': 200}]
        fit = Fit.mpfit_do(g4res, xs, prof, parguess)'''


    fa = {"x": x, "y": y}
    if error is not None:
        fa["error"] = error

    lsf = mpfit.mpfit(residual_fun, parinfo=parinfo, functkw=fa,
            quiet=quiet, maxiter=maxiter)


    return lsf



