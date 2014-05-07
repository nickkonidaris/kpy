# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\npk\.spyder2\.temp.py
"""


import pyzdde.zdde as pyz

def create_genf(link, row, samp=1, wave=0, field=1, dist=13.5,
                typ=1, refp=1, no_diff_lim=1):
    '''Push a GENF onto the merit function at *row*
    
    Inputs:
        link: The PyZDDE object connected to ZEMAX
        
        row [int]: Row number
        
        Remaining parameters: see notes from manual
    
    Returns:
        None
    
    Notes from ZEMAX Manual:
    This operand computes the fraction of geometric encircled, ensquared, 
    x only, or y only (enslitted) energy at a given distance from the 
    reference point defined by Dist. 
    For focal mode, the units are 
    micrometers. For afocal mode, the units are afocal mode units. The other 
    parameters are:
    
    Samp: The pupil sampling, where 1 yields 32 x 32, 2 yields 64 x 64 etc.
        
    Wave: The wavelength number to use (use 0 for polychromatic).
        
    Field: The field number.
        
    Type: 1 for encircled, 2 for x only, 3 for y only, and 4 for ensquared.
        
    Refp: The reference point to use. Use 0 for chief ray, 1 for centroid, 2 for vertex, and 3 for middle of the spot.
            
    No Diff Lim: If 0, the results are scaled by the diffraction limit, otherwise, no accounting of diffraction is done.

    '''
    
    link.zInsertMFO(1)
    link.zSetOperand(row, 1, 'GENF')
    link.zSetOperand(row, 2, samp)
    link.zSetOperand(row, 3, field)
    link.zSetOperand(row, 4, dist)
    link.zSetOperand(row, 5, typ)
    link.zSetOperand(row, 6, refp)
    link.zSetOperand(row, 7, no_diff_lim)

def set_field_angle(link, num, angleX, angleY, weight=1.0):
    '''Set the field angle _num_ to _angleX_ and _angleY_
    
    Inputs:
        link: PyZDDE object connected to zemax
        num: Field position number
        angleX: Angle in X direction (in angular units)
        angleY: Angle in Y direction
        weight: Field weight
    
    Returns:
        None'''
    
    if num == 0: raise Exception("There is no field angle 0")
        
    link.zSetField(
        num, # field number
        angleX, angleY, # Angles x/y
        weight) # Field weight

if __name__ == '__main__':
    l0 = pyz.PyZDDE()
    status = l0.zDDEInit()
    print status
    
    l0.zLoadFile("C:\\Users\\npk\\Dropbox\\zemax_runs\\2014_apr_1_mobie_a")
        
    import numpy as np
    import pylab as pl
    
    l = []
    
    VALUE = 10 # Merit function value
    create_genf(l0, 1, 1, dist=13.5)
    
    XS = np.arange(-0.07, 0.07, 0.1)
    YS = np.arange(-0.06, 0.06, 0.1)
    res = []
    for aX in XS:
        for aY in YS:
            set_field_angle(l0, 1, aX, aY)
            l0.zOptimize(-1, 1)
            genf = l0.zGetOperand(1,10)
            print aX,aY, genf
            res.append(genf)
            
    