# IAPWS / Wagner & Pruß (2002)
# https://pubs.aip.org/aip/jpr/article/31/2/387/241937/The-IAPWS-Formulation-1995-for-the-Thermodynamic

# Feistel & Wagner (2006)
# https://www.sciencedirect.com/science/article/pii/S0016703706020461

import numpy as np
from scipy.interpolate import PchipInterpolator

# Coefficients
_a1 = -7.85951783
_a2 = 1.84408259
_a3 = -11.7866497
_a4 = 22.6807411
_a5 = -15.9618719
_a6 = 1.80122502

_b1 = 1.99274064
_b2 = 1.09965342
_b3 = -0.510839303
_b4 = -1.75493479
_b5 = -45.5170352
_b6 = -6.74694450e5

_c1 = -2.03150240
_c2 = -2.68302940
_c3 = -5.38626492
_c4 = -17.2991605
_c5 = -44.7586581
_c6 = -63.9201063

# Water properties 
T_boil = 373.0
P_boil = 1.0e5    # Pa

T_trip = 273.16
P_trip = 611.655 # Pa

T_crit = 647.096 
P_crit = 22.064 * 1.0e6
rho_crit = 322.0

mmw = 0.018015  # kg/mol

big_number = 1e99

# Functions 

# saturation pressure (gas phase) [Pa]
def psat_liquid(t:float):

    if t < T_trip:
        raise Exception("Evaluation of liquid psat inside solid region")
    
    t_c = 647.096
    p_c = 22.064 * 1.0e6 # Pa

    q = 1 - t/t_c

    rhs = (t_c/t) * (_a1*q + _a2*q**1.5 + _a3*q**3 + _a4*q**3.5 + _a5*q**4 + _a6*q**7.5)
    out = np.exp(rhs) * p_c
    return out 

# saturation pressure (solid phase) [Pa]
def psat_solid(t:float):

    if t > T_trip:
        raise Exception("Evaluation of solid psat inside liquid region")
    
    t = max(t, 20)

    t_n = 273.16 # K
    p_n = 0.000611657 * 1e6  # Pa

    theta = t/t_n 
    rhs = -13.928169*(1-theta**(-1.5)) + 34.7078238*(1-theta**(-1.25))
    out = np.exp(rhs) * p_n
    return out 


# saturation pressure (either phase) [Pa]
def psat_both(t:float):

    if t > T_crit:
        return big_number
    
    if t > T_trip:
        return psat_liquid(t)
    else:
        return psat_solid(t)
    

# derivative of saturation pressure
def dpsat_dt(t:float):
    
    psat = psat_both(t)

    q = 1 - t/T_crit

    paren = np.log(psat/P_crit) + _a1 + 1.5*_a2*q**0.5 + 3*_a3*q**2 + 3.5*_a4*q**2.5 + 4*_a5*q**3 + 7.5*_a6*q**6.5
    out = -1.0 * (psat/t) * paren 
    return out 

# density of liquid
def rho_liq(t):
    q = 1 - t/T_crit

    rhs = 1 + _b1*q**(1/3.0) + _b2*q**(2/3.0) + _b3*q**(5/3.0) + _b4*q**(16/3.0) + _b5*q**(43/3.0) + _b6*q**(110/3.0)
    out = rhs * rho_crit 
    return out 

# density of vapour
def rho_vap(t):
    q = 1 - t/T_crit

    rhs = _c1 * q**(2/6.0) + _c2*q**(4/6.0) + _c3*q**(8/6.0) + _c4*q**(18/6.0) + _c5*q**(37/6) + _c6*q**(71/6)
    out = np.exp(rhs) * rho_crit 
    return out 


# enthalpy of phase change [J/kg]
def delta_vap(t):

    if t >= T_crit:
        return 0.0

    dpdt = dpsat_dt(t)
    h_liq = 1.0e-6 * 1e3 * t / rho_liq(t) 
    h_vap = 1.0e-6 * 1e3 * t / rho_vap(t)

    out = dpdt * (h_vap - h_liq) * 1e3 

    return out 

# sublimation data [K, kJ/kg]
tab_sublim = [  [T_trip , 2834.4],  # set to T_trip for numerics
                [270 , 2835.2],
                [265 , 2836.3],
                [260 , 2837.2],
                [255 , 2837.9],
                [250 , 2838.4],
                [245 , 2838.7],
                [240 , 2838.8],
                [235 , 2838.7],
                [230 , 2838.5],
                [225 , 2838.0],
                [220 , 2837.4],
                [215 , 2836.6],
                [210 , 2835.6],
                [205 , 2834.5],
                [200 , 2833.2],
                [195 , 2831.7],
                [190 , 2830.0],
                [185 , 2828.1],
                [180 , 2826.1],
                [175 , 2823.9],
                [170 , 2821.6],
                [165 , 2819.0],
                [160 , 2816.3],
                [155 , 2813.5],
                [150 , 2810.4],
                [145 , 2807.2],
                [140 , 2803.9],
                [135 , 2800.3],
                [130 , 2796.6],
                [125 , 2792.7],
                [120 , 2788.7],
                [115 , 2784.5],
                [110 , 2780.1],
                [105 , 2775.5],
                [100 , 2770.7],
                [95  , 2765.7],
                [90  , 2760.6],
                [85  , 2755.2],
                [80  , 2749.6],
                [75  , 2743.9],
                [70  , 2737.9],
                [65  , 2731.7],
                [60  , 2725.3],
                [55  , 2718.7],
                [50  , 2712.0],
                [45  , 2705.0],
                [40  , 2698.0],
                [35  , 2690.8],
                [30  , 2683.5],
                [25  , 2676.2],
                [20  , 2668.8],
                [15  , 2661.6],
                [10  , 2654.6],
                [0   , 2654.6]  # same as 10 K for numerics
                ]

tab_sublim = np.array(tab_sublim)[::-1].T

def delta_sub(t:float):
    itp = PchipInterpolator(tab_sublim[0], tab_sublim[1])
    return itp(t) * 1.0e3
    

def delta_both(t:float):
    if t > T_trip:
        return delta_vap(t)
    else:
        return delta_sub(t)
