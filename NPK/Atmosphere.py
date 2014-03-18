
import scipy as sp
import numpy as np
from scipy.interpolate import interp1d

# Palomar Extinction Data from Hayes & Latham 1975
# (Wavelength in Angstroms, Magnitudes per airmass)
palextinct = [
	(3200, 1.058),
	(3250, 0.911),
	(3300, 0.826),
	(3350, 0.757),
	(3390, 0.719),
	(3448, 0.663),
	(3509, 0.617),
	(3571, 0.575),
	(3636, 0.537),
	(3704, 0.500),
	(3862, 0.428),
	(4036, 0.364),
	(4167, 0.325),
	(4255, 0.302),
	(4464, 0.256),
	(4566, 0.238),
	(4785, 0.206),
	(5000, 0.183),
	(5263, 0.164),
	(5556, 0.151),
	(5840, 0.140),
	(6055, 0.133),
	(6435, 0.104),
	(6790, 0.084),
	(7100, 0.071),
	(7550, 0.061),
	(7780, 0.055),
	(8090, 0.051),
	(8370, 0.048),
	(8708, 0.044),
	(9832, 0.036),
	(10255, 0.034),
	(10610, 0.032),
	(10795, 0.032),
	(10870, 0.031)
]

palextinct = np.array(palextinct)
ext = interp1d(palextinct[:, 0], palextinct[:,1], kind='cubic', bounds_error=False)


# From Turnrose PASP 86 (1974)
# Wavelength [A],  10^-18 erg / sec / cm^2 / Angstrom / Arcsecond^2
skyspec = [
        (3180, 4.30),
        (3260, 5.18),
        (3340, 6.13),
        (3420, 4.75),
        (3500, 4.86),
        (3580, 5.29),
        (3660, 7.24),
        (3740, 4.75),
        (3820, 4.43),
        (3900, 3.45),
        (3980, 4.31),
        (4060, 8.58),
        (4140, 6.09),
        (4220, 5.83),
        (4300, 5.39),
        (4380, 11.40),
        (4460, 6.25),
        (4540, 6.38),
        (4620, 6.16),
        (4700, 6.27),
        (4780, 6.14),
        (4860, 6.45),
        (4940, 6.24),
        (5020, 5.60),
        (5100, 5.80),
        (5180, 6.37),
        (5260, 6.26),
        (5340, 6.56),
        (5420, 7.85),
        (5500, 11.00),
        (5580, 25.40),
        (5660, 7.78),
        (5740, 9.70),
        (5760, 9.43),
        (5920, 11.40),
        (6080, 7.89),
        (6240, 13.00),
        (6400, 9.60),
        (6560, 8.36),
        (6720, 6.67),
        (6880, 9.73),
        (7040, 7.11),
        (7200, 9.53),
        (7360, 13.80),
        (7520, 10.70),
        (7680, 13.20),
        (7840, 23.60),
        (8000, 16.60),
        (8160, 5.54),
        (8320, 22.70),
        (8480, 19.30),
        (8640, 20.10),
        (8800, 36.10),
        (8960, 28.30),
        (9120, 8.22),
        (9280, 21.40),
        (9440, 32.40),
        (9600, 15.80),
        (9760, 26.30),
        (9920, 66.00),
        (10080, 68.30),
        (10240, 99.60),
        (10400, 87.10),
        (10560, 25.80),
        (10720, 64.30),
        (10880, 134.00)
]


skyspec = np.array(skyspec)
skyspec[:,1] *= 1e-18 * 3

# See derivation on pg 83 of SED NB 1 (20 July 2011)
moon_phase = np.array([0., 0.08, 0.16, 0.24, 0.32, 0.40, 0.50])
moon_g = np.array([2e-17, 2.1e-17, 2.15e-17, 2.3e-17, 5.3e-17, 1.7e-16, 3.2e-16])
moon_r = np.array([2.3e-17,2.3e-17,2.3e-17,3.3e-17,3.5e-17,8.3e-17,1.3e-16])
moon_i = np.array([2.8e-17,3.0e-17,3.0e-17,3.3e-17,3.8e-17,7.0e-17,9.0e-17])

sky_ls = (4868., 6290., 7706., 10000)


moon_funs = []
for i in xrange(len(moon_phase)):
    gm = moon_g[i]-moon_g[0]
    rm = moon_r[i]-moon_r[0]
    im = moon_i[i]-moon_i[0]
    zm = im

    ff= np.poly1d(np.polyfit(sky_ls, np.array([gm, rm, im, zm]), 2))

    moon_funs.append(ff)


def sky_function(PHASE):
    '''Returns a function of wavelength that returns photon/s/cm^2/ang/as^2
    
    Args:
        PHASE: Moon Phase 
    
    Returns:
        function(wavelength in angstrom) that returns photon/s/cm^2/ang/as^2

    Example:
        f = sky_function(3)

        print("Sky at 5300 Ang is {0} photon/s/cm^2/ang/as^2".{f(5300)})
    '''


    hc = 1.98644521e-8 # erg angstrom
    sky = skyspec[:,1] + moon_funs[PHASE](skyspec[:,0])
    sky /= hc/skyspec[:,0] # #/s/cm^2/ang/as^2
    skyf = interp1d(skyspec[:,0], sky)

    return skyf


