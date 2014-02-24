import numpy as np
import scipy as sp
import IO

reload(IO)

def cube(cube):
    
    shp = cube['Cube'].shape
    im = np.zeros((shp[1]*shp[2], shp[0]))

    for i in xrange(shp[1]):
        for j in xrange(shp[2]):
            im[i*shp[2] + j,:] = cube['Cube'][:,i,j]

    return im


def ds9(cube):
    
    preamble = '''
# Region file format: DS9 version 4.1
global color=blue dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=1 edit=0 move=0 delete=1 include=1 source=1
physical
line(938,1117,1103,1117) # line=0 0
'''


def ds9(segments):
    preamble = '''
# Region file format: DS9 version 4.1
global color=blue dashlist=8 3 width=1 font="helvetica 6 normal" select=1 highlite=1 dash=0 fixed=1 edit=0 move=0 delete=1 include=1 source=1
physical
'''

    for i in xrange(len(segments)):
        segment = segments[i]

        minx = segment['MinX'][0][0][0]
        maxx = segment['MaxX'][0][0][0]
        pp   = segment['Par1'][0].T[0]
        pf   = np.poly1d(pp)
        offset = segment['MeasuredOffset'][0][0][0]


        preamble += "line(%f,%f,%f,%f) # text={%i}\n" % (minx, pf(minx)+offset, 
                                                        maxx, pf(maxx)+offset, 
                                                        i)

    return preamble



def ds92(segments):
    preamble = '''
# Region file format: DS9 version 4.1
global color=blue dashlist=8 3 width=1 font="helvetica 6 normal" select=1 highlite=1 dash=0 fixed=1 edit=0 move=0 delete=1 include=1 source=1
physical
'''

    for i in xrange(len(segments)):
        segment = segments[i]

        X = segment['BlockCenterX'][0].T[0]
        Y = segment['BlockCenterY'][0].T[0]
        pp   = np.polyfit(X,Y,1)
        pf   = np.poly1d(pp)
        minx = X[0]
        maxx = X[-1]
        offset = 0.0

        preamble += "line(%f,%f,%f,%f) # text={%i}\n" % (minx, pf(minx)+offset, 
                                                        maxx, pf(maxx)+offset, 
                                                        i)

    return preamble





