
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

