import numpy as np

# Antoine equation: fit for Psat as a function of temperature

# Most of these values taken from NIST Webbook

gases = {}

gases["NH3"] = {
    "cite": ["Stull 1947", "Brunner 1988"],
    "fit":  [   
                [164.0,    3.18757,    506.713,     -80.78],    # T_min, A, B, C
                [239.6,    4.86886,    1113.928,    -10.409]    # ^
            ],
    "T_crit": 405.4,  
    "T_trip": 195.4,
}

gases["CO2"] = {
    "cite": ["Giauque and Egan, 1937", "Suehiro et al. 1996", "Marsh 1987"],
    "fit":  [   
                [154.26,    6.81228,    1301.679,   -3.494]
            ],
    "T_crit": 304.18, 
    "T_trip": 216.58, 
}

gases["CH4"] = {
    "cite": ["Prydz and Goodwin 1972"],
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



for k in gases.keys():
    gases[k]["cite"].append("P.J. Linstrom and W.G. Mallard, Eds., NIST Chemistry WebBook, NIST Standard Reference Database Number 69")


def evaluate(t,gas_dict):
    fit = gas_dict["fit"]

    A = []
    for f in fit:
        t_min = f[0]
        if t > t_min:
            A = [f[1],f[2],f[3]]

    # zero saturation pressure at low temperatures
    # requires the gas to condense/solidify
    if len(A) == 0:
        return 0.0

    return 10**( A[0] - (A[1]/(t+A[2])) ) * 1e-5  # Pa