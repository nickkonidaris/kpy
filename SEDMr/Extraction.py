
class Extraction():
    seg_id = None 
    ok = None
    xrange = None
    yrange = None
    poly = None
    spec = None
    hg_lines = None

    def __init__(self, seg_id=None, ok=None, xrange=None, 
        yrange=None, poly=None, spec=None,
        hg_lines = None):
        

        self.seg_id = seg_id
        self.ok = ok
        self.xrange = xrange
        self.yrange = yrange
        self.poly = poly
        self.spec = spec
        self.hg_lines = hg_lines
