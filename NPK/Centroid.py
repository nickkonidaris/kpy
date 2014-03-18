import numpy as np

def wm(box):
    '''Computer centroid via weighted mean'''

    total = np.float(box.sum())
    y = np.arange(box.shape[0])
    x = np.arange(box.shape[1])

    xnum = (x*box.sum(0)).sum()
    xcm = xnum/total
    ynum = (y * box.sum(1)).sum()
    ycm = ynum/total

    #print "x: %f +- %f | y: %f +- %f " % (xcm, xcm * np.sqrt(np.abs(1/xnum) + np.abs(1/total)), ycm, ycm * np.sqrt(np.abs(1/ynum) + np.abs(1/total)))

    


    return np.array([xcm, ycm])
