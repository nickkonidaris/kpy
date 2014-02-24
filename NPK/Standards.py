'''Loads standard spectra from Oke 1990 (Aj, 99, 1621).

See: 
    ftp://ftp.eso.org/pub/stecf/standards/okestan/aaareadme.oke

Defines:
    Standards[standard_name]: [Wavelength, flux, flux, bin]
    units: 'Wavelength [A], Flux [erg/s/cm/cm/A 10**16], flux [mJy], bin [A]'
'''
import os
import numpy as np



units = 'Wavelength [A], Flux [erg/s/cm/cm/A 10**16], flux [mJy], bin [A]'

dir = os.path.join(os.path.dirname(__file__), 'standard_stars')

files = os.listdir(dir)

Standards = {}

for file in files:
    if file[0] != 'f': continue

    std_name = file[1:-4]
    dat = np.loadtxt(os.path.join(dir, file))
    Standards[std_name] = dat

