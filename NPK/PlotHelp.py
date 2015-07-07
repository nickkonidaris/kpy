import numpy as np
import pylab as pl


def transparent_legend(legend, alpha=.5):
    
    leg = pl.legend(legend, fancybox=True)
    leg.get_frame().set_alpha(alpha)
