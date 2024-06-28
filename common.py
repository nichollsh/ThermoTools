import matplotlib as mpl

font = {'family' : 'Times',
        'size'   : 12}

mpl.rc('font', **font)

import numpy as np
import glob, os
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

def pretty_name(form):
    pretty = ""
    for f in form:
        if f.isnumeric():
            f = r"$_{%d}$"%int(f)
        pretty += f 
    return pretty

