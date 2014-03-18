import Atmosphere as AA
import numpy as np

reload(AA)

def abmag_to_flambda(AB, lam_ang):
    '''Convert AB Magnitude to erg/s/cm^2/ang

    Arg:
        AB: AB Magnitude
        lam_ang: Wavelength [angstrom]

    Return:
        erg/s/cm^2/Ang'''

    c = 2.9979e18 # Angstrom / s

    if not (2000 < lam_ang < 20000):
        print "LAM MAY NOT BE IN Angstrom"

    # Return erg/s/cm^2/ang
    return 10**(-(AB+48.6)/2.5)*c/lam_ang**2


def go(
    source_ab=None,
    leff_ang=5000,
    R=5,
    fwhm_as=2,
    phase=0,
    t_s=1,
    Atel_cm2=100,
    n_pix=9,
    RN=5,
    eta=.2,
    coadds=1,
    sky_level=1.0):

    '''Calculates the signal to noise of a source with source_ab

    Args:
        source_ab: Source magnitude in AB
        leff_ang: Central wavelength in Angstrom
        R: Spectral resolution of an element
        fwhm_as: Extraction FWHM in as
        phase: Moon phase
        t_s: Exposure time in second
        Atel_cm2: Telescope area in cm2 
        n_pix: Number of pixels that participate
        RN: The read noise in electron
        eta: The efficiency of the instrument
        coadds: The number of coadds [1]

    Returns:
        {All Arg parameters and
        'epp': Energy per photon in erg
        'dlambda': The full with of the band
        's_npp': The number of photons receive from the object
        'sky_area': The area of extraction'
        'k_npp': The number of sky photons received in the extraction area
        'k_npp_pix': The number of sky photons received in a pixel
        'k_npp_as2': The number of sky photons received per as2
        'noise_obj': The noise from the object [e-]
        'noise_sky': The noise from the sky [e-]
        'noise_read': Total read noise
        'noise': Total number of noise photons
        'snr': The delivered signal to noise
        'sky_level': Multiply sky level by this value
    '''

    results={'AB': source_ab, 'leff': leff_ang, 'R': R, 'fwhm_as': fwhm_as,
        'phase': phase, 't': t_s, 'Atel': Atel_cm2, 'n_pix': n_pix,
        'RN': RN, 'eta': eta, 'coadds': coadds, 'sky_level': sky_level}


    hc = 1.98644521e-8 # erg angstrom
    
    epp = hc/leff_ang
    results['epp'] = epp

    dlambda = leff_ang/R
    results['dlambda'] = dlambda

    t_eff = t_s * coadds
    results['t_eff'] = t_eff

    s_flam = abmag_to_flambda(source_ab, leff_ang) # erg/s/cm2/ang
    s_npp = s_flam/epp * t_eff * Atel_cm2 * dlambda * eta
    results['s_npp'] = s_npp

    skyfun = AA.sky_function(phase)

    sky_area = np.pi*(fwhm_as/2.0)**2
    results['sky_area'] = sky_area
    k_flam = skyfun(leff_ang) # photon/s/cm2/ang/as2
    k_npp = k_flam * t_eff * Atel_cm2 * dlambda * sky_area * eta * sky_level

    results['k_npp'] = k_npp
    results['k_npp_pix'] = k_npp/n_pix
    results['k_npp_pas2'] = k_npp/sky_area

    results['noise_obj'] = np.sqrt(s_npp)
    results['noise_sky'] = np.sqrt(k_npp)
    results['noise_read'] = np.sqrt(RN**2*n_pix*coadds)
    results['noise'] = np.sqrt(results['noise_obj']**2 + 
        results['noise_sky']**2 +
        results['noise_read']**2)
    results['snr'] = s_npp/results['noise']

    results['rtel'] = np.sqrt(Atel_cm2/np.pi)
    results['l0'] = leff_ang-dlambda/2
    results['l1'] = leff_ang+dlambda/2

    

    print("                          Signal        Noise")
    print("  Time    Mag   SNR     Star Sky/as2    Star     Sky    Read  R [cm] Waverange  as2")
    print("{t_eff:6.1f} {AB:6.2f} {snr:5.1f}  {s_npp:7.1f} {k_npp_pas2:7.1f} {noise_obj:7.1f} {noise_sky:7.1f} {noise_read:7.1f} {rtel:7.1f} {l0:.0f}-{l1:.0f} {sky_area:4.1f}".format(**results))

    return results



if __name__ == '__main__':

    eta = (.9 * # camera
        .7 *# atmosphere
        .4) #ccd

    Atel_cm2 = np.pi * (45)**2
    go(20.0, t_s=120.0, Atel_cm2=Atel_cm2, eta=0.25, phase=0, leff_ang=5500,
        fwhm_as=1)
    print "Should print out about SNR 33"

    ergpersectowatt = 1e-7
    cm2tom2 = 1e-4
    angtomicron = 1e-4

    toWpm2pmicron = ergpersectowatt / cm2tom2 / angtomicron
    
    # Code Check
    assert( np.abs(abmag_to_flambda(0, 5240) * toWpm2pmicron - 4.069e-8) < .2e-8)
