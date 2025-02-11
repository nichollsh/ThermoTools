import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad
import os

from thermotools import get_gendir
import thermotools.moles as moles

# Antoine equation: fit for Psat as a function of temperature

# Most of these values taken from NIST Webbook

gases = {}

gases["NH3"] = {
    "cite": ["Stull 1947", "Brunner 1988", "Overstreet and Giauque, 1937"],
    "fit":  [   
                [164.0,    3.18757,    506.713,     -80.78],    # T_min, A, B, C
                [239.6,    4.86886,    1113.928,    -10.409]    # ^
            ],
    "T_crit": 405.4,  
    "T_trip": 195.4,
}

gases["CO2"] = {
    "cite": ["Giauque and Egan, 1937", "Suehiro et al. 1996", "Marsh 1987","Stephenson and Malanowski, 1987"],
    "fit":  [   
                [154.26,    6.81228,    1301.679,   -3.494]
            ],
    "T_crit": 304.18, 
    "T_trip": 216.58, 
}

gases["CH4"] = {
    "cite": ["Prydz and Goodwin 1972","Stephenson and Malanowski, 1987"],
    "fit":  [  
                [ 90.99, 3.9895, 443.028, -0.49]
            ],
    "T_crit": 190.6, 
    "T_trip": 90.67,
}

gases["H2"] = {
    "cite": ["van Itterbeek, Verbeke, et al., 1964", "Onnes, Crommelin, et al., 1917", " Clusius and Weigand, 1940"],
    "fit":  [   
                [ 21.01, 3.54314, 99.395, 7.726]
            ],
    "T_crit": 33.18, 
    "T_trip": 13.95,
}

gases["N2"] = {
    "cite": ["Edejer and Thodos, 1967","Jacobsen, Stewart, et al., 1986"],
    "fit":  [  
                [ 63.14 , 	3.7362, 264.651,    -6.788	]
            ],
    "T_crit": 126.19, 
    "T_trip": 63.14	, 
    "h_vap": 6.1,   
}

gases["SO2"] = {
    "cite": ["Stull, 1947","Travers and Usher, 1906","Giauque and Stephenson, 1938"],
    "fit":  [   
                [ 177.7,    3.48586,    668.225,    -72.252],
                [ 263.,	    4.37798,    966.575,    -42.071]
            ],
    "T_crit": 430.34, 
    "T_trip": 197.64, 
}

gases["N2O"] = {
    "cite": ["Stull, 1947","Ohgaki, Umezono, et al., 1990","Giauque and Stephenson, 1938"],
    "fit":  [   
                [ 129.8, 4.37799, 621.077, -44.659	],
            ],
    "T_crit": 309.56, 
    "T_trip": 182.33, 
}

gases["O2"] = {
    "cite": ["Brower and Thodos, 1968","Pentermann and Wagner, 1978	","Henning and Otto, 1936"],
    "fit":  [   
                [ 54.36, 3.9523, 340.024, -4.144],
            ],
    "T_crit": 154.58, 
    "T_trip": 54.33,
}

gases["H2S"] = {
    "cite": ["Stull, 1947","Cubitt, Henderson, et al., 1987","Goodwin, 1983"],
    "fit":  [   
                [ 138.8 , 4.43681,  829.439, -25.412],
                [ 212.8 , 4.52887,  958.587, -0.539]
            ],
    "T_crit": 373.3, 
    "T_trip": 187.66, 
}

gases["O3"] = {
    "cite": ["Stull, 1947","airliquide.com","Jenkins and Birdsall, 1952	"],
    "fit":  [   
                [ 92.8,	4.23637,	712.487,	6.982],
            ],
    "T_crit": 261.15, 
    "T_trip": 80.15, 
}

gases["HCN"] = {
    "cite": ["Sinozaki, Hara, et al., 1926", "Giauque and Ruehrwein, 1939", "https://cameochemicals.noaa.gov/chris/HCN.pdf"],
    "fit":  [   
                [ 256.73,	4.67417,	1340.791,	-11.592],
            ],
    "T_crit": 456.7, 
    "T_trip":259.86, 
}

gases["SiH4"] = {
    "cite": ["Stull, 1947","airliquide.com"],
    "fit":  [   
                [ 93.8 ,	4.22228,	703.987,	5.352],
            ],
    "T_crit": 269.7, 
    "T_trip":88.48, 
}


for k in gases.keys():
    gases[k]["cite"].append("P.J. Linstrom and W.G. Mallard, Eds., NIST Chemistry WebBook, NIST Standard Reference Database Number 69")

def antoine(t:float,gas:str):
    gas_dict = gases[gas]
    fit = gas_dict["fit"]

    A = []
    for f in fit:
        t_min = f[0]
        if t > t_min:
            A = [f[1],f[2],f[3]]

    if len(A) == 0:
        raise Exception("Temperature out of range")

    return 10**( A[0] - (A[1]/(t+A[2])) ) * 1e5  # Pa


# find reference temperature & pressure using Antoine equation
def reference_point(gas:str):
    t_ref = gases[gas]["fit"][0][0] + 20.0
    p_ref = antoine(t_ref, gas)
    return t_ref, p_ref 

# calculate saturation pressure using Clausius-Clapeyron
def cc_psat(t:float, gas:str, mmw:float):

    # find l_vap data
    fname = os.path.join(get_gendir(),"lv","dat")+"/%s.csv"%gas
    t_ref, p_ref = reference_point(gas)
    data = np.loadtxt(fname, delimiter=',').T
    vap_T = data[0]
    vap_L = data[1] * mmw  # convert to J/mol 

    itp = PchipInterpolator(vap_T, vap_L)

    def integrand(t:float):
        L = itp(t)
        return L / (8.314 * t * t)
    
    integ = quad(integrand, t_ref, t)
    rhs = integ[0]
    out = np.exp(rhs) * p_ref 
    return out 

