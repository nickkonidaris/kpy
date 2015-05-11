import numpy as np

'''
Struct:
http://stackoverflow.com/questions/1305532/convert-python-dict-to-object
'''

class Struct(object):
    '''Converts a dictionary into an object

    Example:
        a = {'a': 1, 'b': 2}
        o = Struct(a)
        print o.a, o.b
    
    '''
    def __init__(self, **entries): 
        "s = Struct(***dict)"
        self.__dict__.update(entries)


def floatcompress(data, ndig=14):
    '''Adapted from Finkbeiner IDL routine floatcompress'''

    t = data.dtype
    if not ((t == 'float32') or (t == 'float64')):
        raise Exception("Only works on floating point numbers")

    wzer = np.where(data == 0)
    data[wzer] = 1.0

    log2 = np.ceil(np.log(np.abs(data)) / np.log(2.0))
    mant = np.round(data/2.0**(log2 - ndig))/2.0**ndig
    out = mant*2.0**log2

    out[wzer] = 0.0
    return out
