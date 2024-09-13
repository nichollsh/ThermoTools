import matplotlib as mpl

font = {'family' : 'serif',
        'size'   : 12}

mpl.rc('font', **font)

import matplotlib.pyplot as plt

# latexify a string
def pretty_name(form):
    pretty = ""
    for f in form:
        if f.isnumeric():
            f = r"$_{%d}$"%int(f)
        pretty += f 
    return pretty
