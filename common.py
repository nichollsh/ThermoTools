import matplotlib as mpl

font = {'family' : 'Times',
        'size'   : 11.5}

mpl.rc('font', **font)

import numpy as np
import glob, os
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
