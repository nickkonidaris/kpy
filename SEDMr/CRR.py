
import cosmics
import pyfits as pf
import sys




fn = sys.argv[1]
FF = pf.open(fn)
if FF[0].header['ADCSPEED'] == 2: RN = 21
else: RN = 5

c = cosmics.cosmicsimage(FF[0].data, gain=1, readnoise=RN, sigclip = 5,
    sigfrac=0.5, objlim=4.0)

c.run(maxiter=4, verbose=True)

FF[0].data = c.cleanarray
FF[0].header["CRR"] = ("Processed", "With cosmics.py 0.4")

import os

outn = os.path.join(os.path.dirname(fn), "crr_" + os.path.basename(fn))

FF.writeto(outn)
    
