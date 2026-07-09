import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad
import os

from thermotools import get_gendir, get_inpdir

mmHg = 1/0.0075006156130264  # [Pa] from https://www.convertunits.com/from/mmHg/to/Pa
dyne_per_cm2 = 1e-1  # [Pa] from https://www.convertunits.com/from/dyne/cm2/to/Pa
bar = 1e5  # [Pa] from https://www.convertunits.com/from/bar/to/Pa

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
    "prefer": "cc",
}

gases["CO2"] = {
    "cite": ["Giauque and Egan, 1937", "Suehiro et al. 1996", "Marsh 1987","Stephenson and Malanowski, 1987"],
    "fit":  [   
                [154.26,    6.81228,    1301.679,   -3.494]
            ],
    "T_crit": 304.18, 
    "T_trip": 216.58, 
    "prefer": "antoine",
}

gases["CH4"] = {
    "cite": ["Prydz and Goodwin 1972","Stephenson and Malanowski, 1987"],
    "fit":  [  
                [ 90.99, 3.9895, 443.028, -0.49]
            ],
    "T_crit": 190.6, 
    "T_trip": 90.67, 
    "prefer": "cc",
}

gases["H2"] = {
    "cite": ["van Itterbeek, Verbeke, et al., 1964", "Onnes, Crommelin, et al., 1917", " Clusius and Weigand, 1940"],
    "fit":  [   
                [ 21.01, 3.54314, 99.395, 7.726]
            ],
    "T_crit": 33.18, 
    "T_trip": 13.95,
    "prefer": "cc",
}

gases["N2"] = {
    "cite": ["Edejer and Thodos, 1967","Jacobsen, Stewart, et al., 1986"],
    "fit":  [  
                [ 63.14 , 	3.7362, 264.651,    -6.788	]
            ],
    "T_crit": 126.19, 
    "T_trip": 63.14	, 
    "h_vap": 6.1, 
    "prefer": "cc",  
}

gases["SO2"] = {
    "cite": ["Stull, 1947","Travers and Usher, 1906","Giauque and Stephenson, 1938"],
    "fit":  [   
                [ 177.7,    3.48586,    668.225,    -72.252],
                [ 263.,	    4.37798,    966.575,    -42.071]
            ],
    "T_crit": 430.34, 
    "T_trip": 197.64, 
    "prefer": "cc",
}

gases["N2O"] = {
    "cite": ["Stull, 1947","Ohgaki, Umezono, et al., 1990","Giauque and Stephenson, 1938"],
    "fit":  [   
                [ 129.8, 4.37799, 621.077, -44.659	],
            ],
    "T_crit": 309.56, 
    "T_trip": 182.33, 
    "prefer": "antoine",
}

gases["O2"] = {
    "cite": ["Brower and Thodos, 1968","Pentermann and Wagner, 1978	","Henning and Otto, 1936"],
    "fit":  [   
                [ 54.36, 3.9523, 340.024, -4.144],
            ],
    "T_crit": 154.58, 
    "T_trip": 54.33,
    "prefer": "cc",
}

gases["H2S"] = {
    "cite": ["Stull, 1947","Cubitt, Henderson, et al., 1987","Goodwin, 1983"],
    "fit":  [   
                [ 138.8 , 4.43681,  829.439, -25.412],
                [ 212.8 , 4.52887,  958.587, -0.539]
            ],
    "T_crit": 373.3, 
    "T_trip": 187.66, 
    "prefer": "antoine",
}

gases["O3"] = {
    "cite": ["Stull, 1947","airliquide.com","Jenkins and Birdsall, 1952	"],
    "fit":  [   
                [ 92.8,	4.23637,	712.487,	6.982],
            ],
    "T_crit": 261.15, 
    "T_trip": 80.15, 
    "prefer": "antoine",
}

gases["HCN"] = {
    "cite": ["Sinozaki, Hara, et al., 1926", "Giauque and Ruehrwein, 1939", "https://cameochemicals.noaa.gov/chris/HCN.pdf"],
    "fit":  [   
                [ 256.73,	4.67417,	1340.791,	-11.592],
            ],
    "T_crit": 456.7, 
    "T_trip":259.86, 
    "prefer": "cc",
}

gases["SiH4"] = {
    "cite": ["Stull, 1947","airliquide.com"],
    "fit":  [   
                [ 93.8 ,	4.22228,	703.987,	5.352],
            ],
    "T_crit": 269.7, 
    "T_trip":88.48, 
    "prefer": "antoine",
}

gases["Na"] = {
    "cite": ["Rodebush and Walters, 1930"],
    "fit":  [   
                [ 924.0 ,	2.46077,	1873.728,	-416.372],
            ],
    "T_crit": 2573.0, 
    "T_trip":370.98, 
    "prefer": "antoine",
}

gases["K"] = {
    "cite": ["Fiock and Rodebush, 1926"],
    "fit":  [   
                [679.4, 4.45718, 4691.58, 24.195], 
            ],
    "T_crit": 2223, 
    "T_trip":336.35, 
    "prefer": "antoine",
}

gases["Ca"] = {
    "cite": ["Hartmann and Schneider, 1929"],
    "fit":  [   
                [1254.0,	2.78473,	3121.368,	-594.591], 
            ],
    "T_crit": 5000.0, # not known?
    "T_trip":1112.9, 
    "prefer": "antoine",
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

def woitke_load_gas(gas:str):
    # load gas data from Woitke+2017 derived CSV File
    # these are from Table D2
    fname = os.path.join(get_inpdir(),"psat","web", "woitke2017_TableD2.csv")
    data = np.loadtxt(fname, delimiter=',', dtype=str, comments='#', skiprows=1).T

    # locate row for this gas
    igas = -1
    for i,g in enumerate(data[1]):
        if gas in g:
            print("Querying Woitke+2017 phase change data for %s"%gas)
            if igas>-1:
                print(f"    multiple entries for {gas}, found: {data[1][igas]} and {data[1][i]}")
                if "[s" in data[2][i]:
                    igas = i
            else:
                igas = i

    # not found?
    if igas == -1:
        print("Gas %s not found in Woitke+2017 phase change data (Table D2)"%gas)
        return None
    else:
        print(f"    using {data[1][igas]} data for {gas}")
    
    def _val_or_default(val, default):
        if val == "":
            return default
        else:
            return float(val)

    # convert to dictionary
    data = {
        # metadata
        "gas": str(data[1][igas]).split("[")[0],
        "fit_type": int(data[2][igas]),
        "T_min": _val_or_default(data[4][igas], 0.0),
        "T_max": _val_or_default(data[5][igas], 1e9),

        # density of condensed phase (g/cm^3)
        "density": float(data[3][igas]),

        # fit coefficients
        "c0": _val_or_default(data[6][igas],  0.0),
        "c1": _val_or_default(data[7][igas],  0.0),
        "c2": _val_or_default(data[8][igas],  0.0),
        "c3": _val_or_default(data[9][igas],  0.0),
        "c4": _val_or_default(data[10][igas], 0.0),
    }

    # only store cases with fits for psat, skip deltaG cases
    if data["fit_type"] in [1,2,5]:
        return None

    return data

def woitke_load_all():
    # load all Woitke+2017 data
    fname = os.path.join(get_inpdir(),"psat","web", "woitke2017_TableD2.csv")
    data = np.loadtxt(fname, delimiter=',', dtype=str, comments='#', skiprows=1).T

    parsed_data = {}
    seen_gases = list(gases.keys()) # don't load these twice

    for i in range(len(data[1])):
        gas = str(data[1][i]).split("[")[0]
        if gas not in seen_gases:
            gas_data = woitke_load_gas(gas)
            if gas_data is not None:
                parsed_data[gas] = gas_data
                seen_gases.append(gas)

    return parsed_data

# load the data globally
woitke_data = woitke_load_all()
woitke_gases = list(woitke_data.keys())
print("Loaded Woitke+2017 phase change data for %d gases"%len(woitke_gases))
print("    gases: %s"%(", ".join(woitke_gases)))

def woitke_psat(t:float, gas:str):
    # calculate saturation pressure using Woitke+2017 data

    # get this gas
    if gas not in woitke_gases:
        raise Exception("Gas %s not found in Woitke+2017 phase change data (Table D2)"%gas)
    gas_data = woitke_data[gas]
    c0 = gas_data["c0"]
    c1 = gas_data["c1"]
    c2 = gas_data["c2"]
    c3 = gas_data["c3"]
    c4 = gas_data["c4"]

    # handle temperature range
    t = float(t)
    t = max(t, gas_data["T_min"])
    t = min(t, gas_data["T_max"])

    # celcius for some fits
    tc = t - 273.15

    # these are the various fit equations
    match gas_data["fit_type"]:
        case 1:
            #  Sharp & Huebner (1990)
            raise Exception("Woitke+2017 fit type 1 not yet implemented")
        case 2:
            #  NIST-JANAF data (Chase 1986) 
            raise Exception("Woitke+2017 fit type 2 not yet implemented")
        case 3:
            #  NIST-JANAF data (Chase 1986) 
            ln_P =  c0 / t + c1 + c2 * t + c3 * t**2 + c4 * t**3
            return np.exp(ln_P) * dyne_per_cm2 # dyne/cm^2 to Pa
        case 4:
            #  NIST-JANAF data (Chase 1986) 
            ln_P = c0 + c1 / (t + c2)
            return np.exp(ln_P) * dyne_per_cm2 # dyne/cm^2 to Pa
        case 5:
            #  NIST-JANAF data (Chase 1986) 
            raise Exception("Woitke+2017 fit type 5 not yet implemented")
        case 6:
            #  Yaws (1999)
            log10_P = c0 + c1 / t + c2 * np.log10(t) + c3 * t + c4 * t**2
            return 10**log10_P * mmHg # mmHg to Pa
        case 7:
            #  Ackerman & Marley (2001)
            P = c0 * np.exp( (c1*tc + tc**2/c2) / (tc + c3) )
            return P * dyne_per_cm2 # dyne/cm^2 to Pa
        case 8:
            #  Weast (1971)
            ln_P = c0 + c1 / t + c2 / t**2
            return np.exp(ln_P) * bar # bar to Pa
        case 9:
            #   Prydz & Goodwin (1972)
            log10_P = c0 + c1/(t + c2)
            return 10**log10_P * bar # bar to Pa
        case 10:
            # Gail et al. (2013)
            ln_P = c0/t + c1
            return np.exp(ln_P) * dyne_per_cm2 # dyne/cm^2 to Pa
        
        case _:
            raise Exception("Woitke+2017 fit type %d not yet implemented"%gas_data["fit_type"])
        