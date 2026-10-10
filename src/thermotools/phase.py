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

# Sulfur allotropes over solid and liquid sulfur, from Zahnle et al. (2016) Eqs. 4 and 5.
#     These blend the fits of Lyons (2008) across the melting point, with S6 and S7 lumped into S8.
#     Piecewise fits of the form ln(P/bar) = a - b/T, evaluated by lnfit_psat().
#     The coldest datum is at ~310 K, so lower temperatures are an extrapolation.
gases["S8"] = {
    "cite": ["Zahnle, Marley, Morley, and Moses, 2016, ApJ 824, 137", "Lyons, 2008, J. Sulfur Chem. 29, 269"],
    "fit_ln":  [
                [  0.0,  20.0,  11800.0],   # T_min, a, b
                [413.0,   9.6,   7510.0],   # ^
            ],
    "T_crit": 717.15,   # not given by Zahnle et al.; upper limit of Yaws fit for total sulfur vapour
    "T_trip": 398.0,    # melting point quoted by Zahnle et al.
    "prefer": "lnfit",
}

gases["S2"] = {
    "cite": ["Zahnle, Marley, Morley, and Moses, 2016, ApJ 824, 137", "Lyons, 2008, J. Sulfur Chem. 29, 269"],
    "fit_ln":  [
                [  0.0,  27.0,  18500.0],   # T_min, a, b
                [413.0,  16.1,  14000.0],   # ^
            ],
    "T_crit": 717.15,   # as for S8
    "T_trip": 398.0,    # as for S8
    "prefer": "lnfit",
}

# Antoine equation fits from Yaws & Satyro (2015), Table 2 (inorganic compounds)
# Transcribed from data/lv/web/yaws2015_inorganic.pdf, in the units of the source:
#     log10(P/mmHg) = A - B/(T/degC + C),   valid for TMIN <= T <= TMAX [degC]
# These are converted to the form used by antoine() when added to `gases` (see yaws_load_all).
# Commented-out entries cannot be used by the pipeline: either the formula contains an element
#     without a standard atomic weight in data/mmw/web/elements.txt, or the formula is duplicated.
#   formula:      (A,         B,           C,        TMIN,     TMAX),     #  NO  name
yaws_inorganic = {
#   "Ac":         (  8.85429,  20734.8351,   273.15,  1599.85,  3198.00),  #   1 actinium [element not in mmw table]
    "Ag":         (  8.10019,  12288.3589,   191.93,  1009.85,  2180.00),  #   2 silver
    "AgBr":       (  6.21353,   5179.6010,    52.17,   569.00,  1502.00),  #   3 silver bromide
    "AgCl":       (  8.61748,  10917.5071,   356.11,   912.00,  1547.00),  #   4 silver chloride
    "AgI":        (  6.93625,   6695.6629,   145.03,   820.00,  1506.00),  #   5 silver iodide
    "Al":         (  8.24073,  14524.8870,   192.25,  1209.00,  2517.66),  #   6 aluminum
    "AlB3H12":    (  7.55341,   1419.8359,   259.36,   -52.20,    44.50),  #   7 aluminum borohydride
    "AlBr3":      (  6.73909,   1552.6852,   147.43,    81.30,   255.00),  #   8 aluminum bromide
    "AlCl3":      ( 15.75427,   5649.0504,   258.81,   100.00,   180.00),  #   9 aluminum chloride
    "AlF3":       ( 12.51245,  14983.2501,   279.63,   744.00,  1276.00),  #  10 aluminum fluoride
    "AlI3":       (  6.93225,   2000.0426,   111.66,   178.00,   382.00),  #  11 aluminum iodide
    "Al2O3":      ( 10.96973,  25616.3739,   186.85,  2147.95,  2980.00),  #  12 aluminum oxide
    "Al2S3O12":   ( -0.14253,    651.9161,  -975.34,   572.00,   770.05),  #  13 aluminum sulfate
#   "Am":         (  9.46480,  14359.6923,   273.15,   965.85,  1082.85),  #  14 americium [element not in mmw table]
    "Ar":         (  6.84064,    340.2707,   271.80,  -189.37,  -122.29),  #  15 argon
    "As":         ( 10.60233,   6313.8843,   216.09,   280.00,   601.00),  #  16 arsenic
    "AsBr3":      (  7.22542,   1953.4063,   228.62,    41.80,   221.00),  #  17 arsenic bromide
    "AsCl3":      (  7.44841,   1670.0897,   235.64,   -11.40,   130.00),  #  18 arsenic trichloride
    "AsF3":       (  7.30253,   1189.5299,   211.22,    -2.50,    57.80),  #  19 arsenic trifluoride
    "AsF5":       (  5.74019,    366.9598,   181.14,  -117.90,   -52.80),  #  20 arsenic pentafluoride
    "AsH3":       (  6.88872,    790.0472,   259.60,  -116.92,    99.85),  #  21 arsine
    "AsI3":       (  7.30188,   2551.7378,   209.53,   164.00,   329.50),  #  22 arsenic triiodide
    "As2O3":      (  4.99289,   1001.4246,     6.99,   133.70,   467.15),  #  23 arsenic trioxide
#   "At":         (  9.64649,   3910.8049,   244.23,    88.00,   259.00),  #  24 astatine [element not in mmw table]
    "Au":         (  8.16676,  16053.2805,   180.98,   952.85,  2856.00),  #  25 gold
    "B":          (  9.44507,  26209.9893,   190.34,  2075.00,  3799.00),  #  26 boron
    "BBr3":       (  6.58372,   1114.5997,   210.01,   -45.00,    91.00),  #  27 boron tribromide
    "BCl3":       (  6.82139,    972.2015,   234.22,  -107.00,   178.80),  #  28 boron trichloride
    "BF3":        (  7.29112,    618.3678,   240.01,  -128.37,   -12.25),  #  29 boron trifluoride
    "BH2CO":      (  7.47035,    893.9306,   258.95,  -139.20,   -64.00),  #  30 borine carbonyl
    "BI3":        (  6.93670,   1698.2347,   210.58,    43.00,   200.00),  #  31 boron triiodide
#   "B2D6":       (  6.30635,    512.0855,   241.99,  -154.90,   -92.50),  #  32 diborane, deutero [element not in mmw table]
    "B2F4":       (  7.59073,   1126.3767,   273.15,   -40.00,   -14.00),  #  33 diborane, tetrafluoro
    "B2H5Br":     (  7.79587,   1454.1004,   279.55,   -93.30,    16.30),  #  34 diborane hydrobromide
    "B2H6":       (  6.83515,    654.3293,   257.97,  -165.50,    16.65),  #  35 diborane
    "B3H6N3":     (  7.19573,   1254.5863,   237.76,   -63.00,    53.00),  #  36 borazine
    "B4H10":      (  7.23266,   1190.0772,   255.46,   -90.90,    18.00),  #  37 tetraborane
    "B5H9":       (  6.46328,    948.6930,   204.82,   -40.40,   295.30),  #  38 pentaborane
    "B5H11":      (  9.78620,   2656.8429,   321.75,   -50.20,    63.00),  #  39 pentaborane
    "B10H14":     (  4.60815,    512.9508,    48.92,    60.00,   163.80),  #  40 decaborane
    "Ba":         (  6.47688,   7779.9111,   266.45,   638.00,  1897.00),  #  41 barium
    "BaF2":       (  8.20759,  12116.9768,    20.38,  1368.00,  2250.00),  #  42 barium fluoride
    "BaO":        ( 12.51502,  21900.0000,   273.15,  1950.00,  2000.00),  #  43 barium oxide
    "Be":         (  8.08396,  13696.6102,   148.52,   823.85,  2483.85),  #  44 beryllium
    "BeB2H8":     ( 15.84584,   6333.8132,   398.57,     1.00,    90.00),  #  45 beryllium borohydride
    "BeBr2":      ( 11.58765,   6557.5715,   274.67,   289.00,   474.00),  #  46 beryllium bromide
    "BeCl2":      ( 11.21338,   6220.1613,   264.49,   291.00,   482.00),  #  47 beryllium chloride
    "BeF2":       (  9.03258,   7534.4285,    55.76,   686.00,  1169.00),  #  48 beryllium fluoride
    "BeI2":       ( 10.16571,   5302.6269,   240.89,   283.00,   487.00),  #  49 beryllium iodide
    "Bi":         (  7.90580,   9022.0737,   231.44,   668.00,  1564.00),  #  50 bismuth
    "BiBr3":      (  8.73564,   4109.4999,   248.90,   261.00,   453.00),  #  51 bismuth tribromide
    "BiCl3":      (  7.31869,   2847.8317,   194.71,   230.50,   447.00),  #  52 bismuth trichloride
    "BrFO3":      (  7.43651,   1195.8000,   260.10,   -75.15,    28.85),  #  53 perbromyl fluoride
    "BrF3":       (  7.72974,   1673.9500,   219.48,    25.85,   157.85),  #  54 bromine trifluoride
    "BrF5":       (  7.85138,   1488.6143,   258.99,   -69.30,    40.50),  #  55 bromine pentafluoride
    "Br2":        (  7.27906,   1321.8955,   241.80,    -7.30,   311.00),  #  56 bromine
    "C":          ( 13.52595,  41703.0146,   280.60,  2565.85,  3634.85),  #  57 graphite
    "CCl2O":      (  7.14063,   1063.6278,   242.13,  -127.78,   181.85),  #  58 carbonyl chloride
    "CF2O":       (  7.41143,    765.5881,   252.08,  -111.26,    23.85),  #  59 carbonyl fluoride
    "CH4N2O":     ( 10.45495,   4591.7946,   273.76,    67.50,    94.90),  #  60 urea
    "CH4N2S":     (  7.57639,   2173.4632,   198.99,   181.00,   580.85),  #  61 thiourea
    "CNBr":       (  6.36918,    722.7494,   142.72,    -0.14,    39.94),  #  62 cyanogen bromide
    "CNCl":       (  7.55525,   1251.9891,   254.99,    -6.50,   175.85),  #  63 cyanogen chloride
    "CNF":        (  6.75643,    689.9197,   224.02,   -76.40,   -46.00),  #  64 cyanogen fluoride
    "CO":         (  6.72527,    295.2279,   268.24,  -205.00,  -140.23),  #  65 carbon monoxide
    "OCS":        (  7.00667,    833.0673,   252.06,  -138.80,   105.65),  #  66 carbonyl sulfide, 'COS' in Yaws
    "COSe":       (  7.17236,   1024.3828,   260.38,  -120.00,   -22.00),  #  67 carbon oxyselenide
    "CO2":        (  7.58828,    861.8205,   271.88,   -56.57,    31.04),  #  68 carbon dioxide
    "CS2":        (  7.21790,   1303.7883,   254.39,  -111.57,   278.85),  #  69 carbon disulfide
    "C2N2":       (  7.04339,    770.7952,   206.27,   -95.80,   -21.10),  #  70 cyanogen
    "C3S2":       (  7.51838,   1986.8195,   250.50,    14.00,   130.80),  #  71 carbon subsulfide
    "Ca":         (  7.10583,   6959.2899,   162.92,   591.00,  1482.00),  #  72 calcium
    "CaCl2":      (  9.22000,  13570.0000,   273.15,   836.85,  1007.85),  #  73 calcium chloride
    "CaF2":       (  7.47685,  10678.0076,  -210.04,  1417.85,  2533.35),  #  74 calcium fluoride
    "Cd":         (  7.66023,   4757.6158,   228.44,   119.85,   767.00),  #  75 cadmium
    "CdBr2":      ( 10.71437,   9463.5051,   364.07,   373.00,   844.00),  #  76 cadmium bromide
    "CdCl2":      (  9.79511,   9774.8935,   453.72,   618.00,   960.00),  #  77 cadmium chloride
    "CdF2":       (  8.15105,   9504.5987,    55.45,  1112.00,  1748.00),  #  78 cadmium fluoride
    "CdI2":       (  7.72916,   4256.7258,   135.98,   296.00,   742.00),  #  79 cadmium iodide
    "CdO":        (  9.13973,  11107.3677,   215.65,  1000.00,  1559.00),  #  80 cadmium oxide
    "CdS":        ( 10.87780,  10021.4278,   273.15,   930.00,   980.00),  #  81 cadmium sulfide
    "Ce":         (  8.54474,  20813.1074,   231.68,  1719.00,  3443.00),  #  82 cerium
    "CeBr3":      ( 11.23239,  14449.4770,   273.15,  1407.00,  1457.00),  #  83 cerium bromide
    "CeCl3":      ( 10.38288,  15005.2628,   273.15,  1677.00,  1727.00),  #  84 cerium chloride
    "ClF":        (  6.90617,    644.1798,   249.77,  -158.15,    26.85),  #  85 chlorine monofluoride
    "ClFO2":      (  6.67696,    809.7800,   218.96,   -82.15,    21.85),  #  86 chloryl fluoride
    "ClFO3":      (  7.08309,    854.6155,   250.03,  -147.74,    95.25),  #  87 perchloryl fluoride
    "ClF3":       (  7.34481,   1086.4081,   231.62,   -73.70,    11.75),  #  88 chlorine trifluoride
    "ClF5":       (  6.87794,    896.9883,   237.51,   -79.20,    24.80),  #  89 chlorine pentafluoride
    "ClNO2":      (  5.37230,    395.4000,   174.00,   -87.15,   -21.15),  #  90 nitryl chloride
    "ClO2":       (  7.29962,   1033.2300,   222.93,   -59.60,   191.85),  #  91 chlorine dioxide
    "Cl2":        (  7.06306,    906.7031,   250.83,  -101.03,   144.00),  #  92 chlorine
    "Cl2O":       (  7.70028,   1296.5974,   266.83,   -98.50,     2.20),  #  93 chlorine monoxide
    "Cl2O4":      (  7.53867,   1404.1800,   257.00,   -44.15,    74.85),  #  94 chlorine perchlorate
    "Cl2O6":      (  5.83238,   1400.9600,   273.15,   151.50,   201.50),  #  95 chlorine hexoxide
    "Cl2O7":      (  7.59797,   1583.6401,   253.72,   -45.30,    82.00),  #  96 chlorine heptoxide
    "Co":         (  8.85240,  18477.2528,   166.22,  1517.00,  2925.00),  #  97 cobalt
    "CoCl2":      (  7.04947,   4946.1939,   137.52,   770.00,  1049.00),  #  98 cobalt chloride
    "CoNC3O4":    (  3.05780,      9.5159,     6.92,    -1.30,    80.00),  #  99 cobalt nitrosyl tricarbonyl
    "Cr":         (  8.48706,  15307.9218,    59.51,  1383.00,  2671.00),  # 100 chromium
    "CrC6O6":     ( 10.71935,   3410.1719,   280.20,    36.00,   151.00),  # 101 chromium carbonyl
    "CrCl2":      ( 12.16568,  12935.2174,   273.15,  1070.00,  1120.00),  # 102 chromium chloride
    "CrCl3":      ( 10.96666,  12720.2532,   273.15,  1250.00,  1300.00),  # 103 chromium chloride
    "CrO2Cl2":    (  6.76400,   1231.6008,   200.16,   -18.40,   117.00),  # 104 chromium oxychloride
    "CrO3":       ( 20.14000,  10300.0000,   273.15,   174.85,   194.85),  # 105 chromium oxide
    "Cs":         (  6.67520,   3515.7582,   255.57,    21.85,   671.00),  # 106 cesium
    "CsBr":       (  7.70150,   7132.1044,   179.48,   748.00,  1300.00),  # 107 cesium bromide
    "CsCl":       (  7.74073,   7233.7811,   191.46,   744.00,  1297.00),  # 108 cesium chloride
    "CsF":        (  7.54319,   6579.0903,   160.10,   712.00,  1251.00),  # 109 cesium fluoride
    "CsI":        (  7.64311,   6849.9814,   158.38,   738.00,  1280.00),  # 110 cesium iodide
    "Cu":         (  8.30991,  15006.7869,   202.14,  1216.00,  2562.00),  # 111 copper
    "CuBr":       (  5.01163,   2858.6991,    -3.41,   572.00,  1345.00),  # 112 cuprous bromide
    "CuCl":       (  4.45573,   2327.4421,   -12.18,   429.85,  1490.00),  # 113 cuprous chloride
    "CuCl2":      ( -4.76869,   3645.4156,  -999.92,   309.70,   521.00),  # 114 cupric chloride
    "CuI":        (  4.81291,   2646.3848,    33.70,   610.00,  1336.00),  # 115 copper iodide
#   "DBr":        (  6.16238,    505.6800,   220.60,   -88.15,   -51.15),  # 116 deuterium bromide [element not in mmw table]
#   "DCN":        (  6.25464,    689.1316,   178.31,   -68.90,    26.20),  # 117 deuterium cyanide [element not in mmw table]
#   "DCl":        (  6.93596,    668.2000,   249.50,  -115.15,   -62.15),  # 118 deuterium chloride [element not in mmw table]
#   "DF":         (  7.21704,   1268.3700,   273.87,   -84.15,    49.85),  # 119 deuterium fluoride [element not in mmw table]
#   "DI":         (  5.60474,    414.6800,   187.87,   -55.15,    -7.15),  # 120 deuterium iodide [element not in mmw table]
#   "D2":         (  6.31539,     92.2133,   276.35,  -254.42,  -234.80),  # 121 deuterium [element not in mmw table]
#   "D2O":        (  8.06347,   1708.9717,   228.35,     3.82,   370.74),  # 122 deuterium oxide [element not in mmw table]
    "Dy":         (  6.67690,   9758.5876,     3.70,  1105.00,  2567.00),  # 123 dysprosium
    "Er":         (  7.07492,  12618.2924,   140.57,  1231.00,  2868.00),  # 124 erbium
    "Eu":         (  6.80515,   6548.9117,   143.35,   590.00,  1523.00),  # 125 europium
    "FHO3S":      (  7.39950,   1521.0000,   174.00,    64.85,   192.85),  # 126 fluorosulfonic acid
    "FNO":        (  7.88277,    986.5042,   257.12,  -132.00,   -59.90),  # 127 nitrosyl fluoride
    "FNO2":       (  6.94723,    685.0522,   240.87,  -156.00,   -72.40),  # 128 nitryl fluoride
    "FNO3":       (  6.53787,    721.4837,   243.29,  -160.00,   -46.00),  # 129 fluorine nitrate
    "F2":         (  6.89461,    320.1001,   267.95,  -219.67,  -128.84),  # 130 fluorine
    "F2O":        (  7.29316,    580.8480,   276.45,  -198.15,   -63.15),  # 131 fluorine oxide
    "F2O2":       (  6.77902,    756.3900,   250.16,  -125.15,   -80.15),  # 132 difluorine dioxide
    "Fe":         (  8.37930,  16213.5121,    88.52,  1435.00,  2860.20),  # 133 iron
    "FeC5O5":     (  9.24491,   2718.2717,   324.13,    -6.50,   103.00),  # 134 iron pentacarbonyl
    "FeCl2":      (  7.23943,   4667.3502,    46.98,   700.00,  1023.85),  # 135 ferrous chloride
    "FeCl3":      ( 28.93093,  32484.7863,   928.01,   108.00,   319.00),  # 136 ferric chloride
    "FeI2":       (  7.09897,   4398.3915,   107.58,   587.00,   926.85),  # 137 iron iodide
#   "Fr":         (  6.28464,   3099.3568,   237.55,   121.00,   673.00),  # 138 francium [element not in mmw table]
    "Ga":         (  8.24091,  13438.1916,   259.39,  1037.00,  2245.00),  # 139 gallium
    "GaBr3":      (  7.64043,   2477.2199,   205.34,   121.50,   226.85),  # 140 gallium bromide
    "GaCl3":      (  7.54055,   2005.3675,   229.36,    77.75,   420.85),  # 141 gallium trichloride
    "Gd":         (  8.09461,  18125.4084,   210.60,  1563.00,  3262.00),  # 142 gadolinium
    "Ge":         (  8.39437,  16953.3407,   240.64,  1371.00,  2833.00),  # 143 germanium
    "GeBr4":      (  6.92749,   1650.1032,   221.42,    43.30,   186.35),  # 144 germanium bromide
    "GeCl4":      (  6.27456,    971.2957,   199.65,   -45.00,    86.55),  # 145 germanium chloride
    "GeHCl3":     (  7.29569,   1304.6359,   220.21,   -41.30,    75.30),  # 146 germane trichloro
    "GeH4":       (  7.08333,    732.4097,   262.43,  -165.89,    34.85),  # 147 germane
    "GeS":        (  6.81801,   6765.5528,   214.38,   615.00,   726.85),  # 148 germanium sulfide
    "Ge2H6":      (  6.98641,   1189.4193,   258.91,   -88.70,    30.80),  # 149 digermane
    "Ge3H8":      (  7.19898,   1588.2499,   257.31,   -36.90,   110.50),  # 150 trigermane
    "HBr":        (  6.93857,    704.4329,   240.30,   -88.00,    90.00),  # 151 hydrogen bromide
    "HCN":        (  8.08422,   1655.8693,   292.23,   -13.32,   183.50),  # 152 hydrogen cyanide
    "HCl":        (  7.25303,    764.8366,   259.93,  -114.18,    51.50),  # 153 hydrogen chloride
    "HClO3S":     (  7.47025,   1749.2078,   227.29,   -80.00,   426.85),  # 154 chlorosulfonic acid
    "HClO4":      (  6.71471,   1348.2463,   239.81,  -101.20,   357.85),  # 155 perchloric acid
    "HF":         (  8.12392,   1680.6466,   301.02,   -83.36,   188.00),  # 156 hydrogen fluoride
    "HI":         (  6.99067,    881.4677,   250.08,   -50.77,   150.70),  # 157 hydrogen iodide
    "HNO3":       (  7.56453,   1431.7460,   222.69,   -41.60,   102.95),  # 158 nitric acid
    "HN3":        (  6.65340,    979.1881,   223.46,   -79.00,    35.70),  # 159 hydrazoic acid
    "HPO2F2":     (  6.73530,   1342.9000,   232.00,     4.85,   151.85),  # 160 difluorophosphoric acid
    "H2":         (  6.14858,     80.9480,   277.53,  -259.20,  -239.97),  # 161 hydrogen
    "H2O":        (  8.05573,   1723.6425,   233.08,     0.01,   373.98),  # 162 water
    "H2O2":       (  8.07489,   1958.1824,   226.80,    -0.41,   457.00),  # 163 hydrogen peroxide
    "H2S":        (  7.11958,    802.2266,   249.61,   -85.47,   100.38),  # 164 hydrogen sulfide
    "H2SO4":      (  8.28538,   2976.8370,   213.95,    25.00,   336.85),  # 165 sulfuric acid
    "H2S2":       (  7.02654,   1242.2189,   228.94,   -40.00,    70.70),  # 166 hydrogen disulfide
    "H2S3":       (  6.80710,   1488.0000,   209.00,    54.85,   200.85),  # 167 dihydrogen trisulfide
    "H2S4":       (  6.94510,   1772.0000,   196.00,   110.85,   273.85),  # 168 dihydrogen tetrasulfide
    "H2S5":       (  7.32010,   2104.0000,   189.00,   152.85,   318.85),  # 169 dihydrogen pentasulfide
    "H2Se":       (  6.46297,    587.4538,   205.99,  -115.30,   -42.00),  # 170 hydrogen selenide
    "H2Te":       (  5.90921,    586.0450,   195.52,   -96.40,    -2.00),  # 171 hydrogen telluride
    "H3NO":       ( 11.71366,   3155.6880,   247.42,    33.10,   109.85),  # 172 hydroxylamine
    "H3SNO3":     (  0.11312,    684.4925,   273.15,    20.00,   100.00),  # 173 sulfamic acid
#   "He":         (  4.82950,      8.5576,   274.36,  -272.15,  -269.83),  # 174 helium-3 [duplicate: helium-3; helium-4 used instead]
    "He":         (  5.27120,     13.5171,   274.58,  -271.39,  -267.95),  # 175 helium-4
    "Hf":         (  9.06491,  30232.9115,   285.82,  2416.00,  4603.00),  # 176 hafnium
    "HfCl4":      ( 11.71000,   5197.0000,   273.15,   202.85,   310.00),  # 177 hafnium chloride
    "HfI3":       ( 19.56000,  10700.0000,   273.15,   301.85,   323.85),  # 178 hafnium iodide
    "Hg":         (  7.89458,   3147.6142,   271.10,   -38.84,  1461.85),  # 179 mercury
    "HgBr2":      (  5.56804,    974.6103,    40.68,   136.50,   322.00),  # 180 mercuric bromide
    "HgCl2":      ( 12.19380,   6369.7914,   379.97,   136.20,   304.00),  # 181 mercuric chloride
    "HgI2":       (  7.11921,   2015.6871,   121.58,   157.50,   354.00),  # 182 mercuric iodide
    "Ho":         (  6.66154,  10226.4211,     4.88,  1159.00,  2700.00),  # 183 holmium
    "ICl":        (  7.69746,   1515.4000,   216.80,    13.85,   123.85),  # 184 iodine chloride
    "IF5":        (  7.24100,   1387.0000,   213.62,     9.85,   134.85),  # 185 iodine pentafluoride
    "IF7":        (  8.58999,   1562.5974,   268.90,   -87.00,     4.80),  # 186 iodine heptafluoride
    "I2":         (  7.10736,   1556.9694,   183.97,   -31.15,   546.00),  # 187 iodine
    "In":         (  7.98843,  11856.7391,   249.38,   900.00,  2072.00),  # 188 indium
    "InBr":       ( 19.02460,  15000.0000,   273.15,   626.00,   656.00),  # 189 indium bromide
    "InBr3":      (  5.89173,    727.4250,  -183.70,   304.60,   364.80),  # 190 indium bromide
    "InCl":       (  8.03000,   4640.0000,   273.15,   211.00,   576.85),  # 191 indium chloride
    "InCl3":      ( 13.62000,   8270.0000,   273.15,   226.85,   426.85),  # 192 indium chloride
    "Ir":         (  9.63425,  30934.6262,   190.68,  2440.00,  4386.00),  # 193 iridium
    "IrF6":       (  8.56221,   1717.5557,   248.71,   -88.00,    53.60),  # 194 iridium fluoride
    "K":          (  6.97431,   4143.9271,   254.47,    63.20,  1949.85),  # 195 potassium
    "KBr":        (  7.60549,   7291.2713,   163.23,   795.00,  1380.00),  # 196 potassium bromide
    "KCl":        (  6.93314,   5810.3508,    17.98,   770.85,  3196.85),  # 197 potassium chloride
    "KF":         (  7.61993,   7726.2515,   128.32,   885.00,  1502.00),  # 198 potassium fluoride
    "KI":         (  6.57543,   5053.2934,    22.90,   745.00,  1344.85),  # 199 potassium iodide
    "KOH":        (  7.08223,   6300.4849,   172.76,   405.85,  1326.85),  # 200 potassium hydroxide
    "Kr":         (  7.01408,    511.6161,   277.13,  -157.37,   -63.80),  # 201 krypton
    "La":         (  8.57782,  21096.9943,   239.17,  1650.00,  3464.00),  # 202 lanthanum
    "LaBr3":      ( 10.98826,  15000.0000,   273.15,  1547.00,  1577.00),  # 203 lanthanum bromide
    "LaCl3":      ( 14.66261,  15000.0000,   273.15,   970.00,  1000.00),  # 204 lanthanum chloride
    "Li":         (  7.80332,   7953.8982,   276.84,   524.30,  1337.10),  # 205 lithium
    "LiBr":       (  7.57994,   6947.3452,   168.43,   748.00,  1310.00),  # 206 lithium bromide
    "LiCl":       (  7.53439,   7299.6503,   185.61,   783.00,  1383.00),  # 207 lithium chloride
    "LiF":        (  8.71722,  11054.9564,   221.14,  1047.00,  1673.00),  # 208 lithium fluoride
    "LiI":        (  8.49354,   7459.8616,   155.25,   723.00,  1173.85),  # 209 lithium iodide
    "Lu":         (  7.83211,  17422.3777,   116.75,  1600.00,  3402.00),  # 210 lutetium
    "Mg":         (  7.16078,   5253.0612,   137.71,   428.00,  1088.00),  # 211 magnesium
    "MgBr2":      (  7.38742,   6095.1378,   104.02,   711.00,  1226.85),  # 212 magnesium bromide
    "MgCl2":      (  7.13801,   6696.7751,   161.05,   778.00,  1412.00),  # 213 magnesium chloride
    "MgF2":       (  7.87219,  12036.3932,    82.67,  1262.85,  2126.85),  # 214 magnesium fluoride
    "MgI2":       (  6.92855,   4121.1468,   109.39,   634.00,   876.85),  # 215 magnesium iodide
    "MgO":        (  6.92276,  13862.2156,  -170.46,  2831.85,  5676.85),  # 216 magnesium oxide
    "Mn":         (  7.50506,   9853.7049,    68.23,   955.00,  2060.00),  # 217 manganese
    "MnCl2":      (  7.03063,   5467.6217,   127.56,   736.00,  1190.00),  # 218 manganese chloride
    "MnF2":       ( 16.16514,  15000.0000,   273.15,   826.00,   856.00),  # 219 manganese fluoride
    "Mo":         (  8.64938,  26943.1312,    31.68,  2400.00,  4950.00),  # 220 molybdenum
    "MoCl5":      (  6.29946,   1850.0000,   273.15,   218.00,   268.00),  # 221 molybdenum chloride
    "MoF5":       (  8.49628,   2709.7945,   268.96,    86.60,   213.60),  # 222 molybdenum fluoride
    "MoF6":       (  7.06989,   1018.5440,   209.14,   -65.50,    34.00),  # 223 molybdenum fluoride
    "MoOF4":      (  8.53528,   2492.6355,   254.83,   -21.00,   186.00),  # 224 molybdenum oxytetrafluoride
    "MoO3":       (  6.05246,   2898.1602,  -233.00,   801.00,  1140.00),  # 225 molybdenum oxide
    "NClO":       (  7.83560,   1314.8605,   270.75,   -59.60,   167.50),  # 226 nitrosyl chloride
    "NCl3":       (  7.00316,   1212.7111,   223.18,   -27.00,    94.00),  # 227 nitrogen trichloride
#   "ND3":        (  6.83776,    746.0100,   221.53,   -74.00,   -33.00),  # 228 heavy ammonia [element not in mmw table]
    "NF3":        (  7.01648,    559.5184,   264.35,  -206.79,   -39.30),  # 229 nitrogen trifluoride
    "NF3S":       (  6.88510,    888.5000,   245.00,   -94.15,     1.85),  # 230 nitrogen fluoride sulfide
    "NH3":        (  7.58743,   1013.7815,   248.83,   -77.74,   132.50),  # 231 ammonia
    "NH4Br":      (  9.74168,   4588.8568,   272.85,   198.30,   396.00),  # 232 ammonium bromide
    "NH4Cl":      (  8.93371,   3215.2979,   192.35,   280.00,   520.00),  # 233 ammonium chloride
    "NH4I":       ( 10.40721,   5286.5838,   297.41,   210.90,   405.00),  # 234 ammonium iodide
    "NH5O":       (  7.47687,   1118.7172,   245.15,   -70.00,    80.00),  # 235 ammonium hydroxide
    "NH5S":       (  8.96559,   1598.3103,   229.34,   -51.10,    33.30),  # 236 ammonium hydrogensulfide
    "NO":         (  7.70299,    473.0296,   249.87,  -163.65,   -93.00),  # 237 nitric oxide
    "NO2":        (  9.74325,   2289.6342,   312.65,   -11.20,   158.20),  # 238 nitrogen dioxide
    "N2":         (  6.72531,    285.5727,   270.09,  -210.00,  -147.05),  # 239 nitrogen
    "N2F4":       (  7.41740,    952.4783,   284.15,  -161.50,    36.20),  # 240 tetrafluorohydrazine
    "N2H4":       (  7.82877,   1698.5808,   229.79,     1.53,   380.00),  # 241 hydrazine
    "N2H4C":      (  9.16524,   1644.6341,   230.00,   -50.60,    31.70),  # 242 ammonium cyanide
    "N2H6CO2":    ( 10.32341,   2267.7036,   245.93,   -26.10,    58.30),  # 243 ammonium carbamate
    "N2O":        (  6.81488,    550.6045,   228.44,   -90.85,    36.42),  # 244 nitrous oxide
    "N2O3":       (  7.60164,   1109.3200,   232.98,  -103.15,   151.85),  # 245 nitrogen trioxide
    "N2O4":       (  4.09407,     97.8287,    59.47,   -11.25,    47.50),  # 246 nitrogen tetraoxide
    "N2O5":       ( 11.93250,   2617.6543,   256.19,   -36.80,    33.00),  # 247 nitrogen pentoxide
    "Na":         (  7.33576,   5060.9675,   254.35,   280.60,   926.85),  # 248 sodium
    "NaBr":       (  7.54103,   7142.9394,   142.08,   746.85,  1446.85),  # 249 sodium bromide
    "NaCN":       (  7.38628,   7836.8557,   243.41,   563.70,  1496.00),  # 250 sodium cyanide
    "NaCl":       (  7.96318,   8432.2350,   194.07,   800.75,  1465.05),  # 251 sodium chloride
    "NaF":        (  8.61974,  11408.3735,   278.32,   995.85,  1786.85),  # 252 sodium fluoride
    "NaI":        (  7.89004,   7361.9912,   165.84,   767.00,  1303.85),  # 253 sodium iodide
    "NaOH":       (  7.46511,   7606.4353,   280.16,   500.00,  1600.00),  # 254 sodium hydroxide
    "Na2SO4":     ( 11.40104,  15657.1880,   273.15,  1514.50,  1564.50),  # 255 sodium sulfate
    "Nb":         (  9.92532,  35189.3422,   251.22,  2669.00,  4740.00),  # 256 niobium
    "NbBr5":      (  9.33000,   4085.0000,   273.15,   254.00,   360.00),  # 257 niobium bromide
    "NbCl4":      ( 12.30000,   6870.0000,   273.15,   303.85,   377.85),  # 258 niobium chloride
    "NbCl5":      (  8.37000,   2870.0000,   273.15,   204.70,   246.85),  # 259 niobium chloride
    "NbF5":       (  6.88947,   1791.5982,   217.93,    86.30,   229.00),  # 260 niobium fluoride
    "Nd":         (  6.95222,  12884.3059,    97.12,  1322.30,  3063.00),  # 261 neodymium
    "NdBr3":      ( 11.15371,  15000.0000,   273.15,  1510.00,  1540.00),  # 262 neodymium tribromide
    "NdCl3":      ( 10.88871,  15000.0000,   273.15,  1570.00,  1600.00),  # 263 neodymium trichloride
    "Ne":         (  6.74431,    115.3285,   275.91,  -248.59,  -228.75),  # 264 neon
    "Ni":         (  8.74585,  17882.3827,   134.99,  1510.00,  2911.00),  # 265 nickel
    "NiC4O4":     (  7.64987,   1479.9228,   267.82,   -23.00,    42.50),  # 266 nickel carbonyl
    "NiCl2":      ( 10.25227,   8220.2197,   130.14,   534.00,   985.00),  # 267 nickel chloride
    "NiF2":       (  8.31881,   7276.8348,  -160.62,  1077.00,  1283.00),  # 268 nickel fluoride
#   "Np":         (  7.90387,  22003.2016,   273.15,  1726.85,  2226.85),  # 269 neptunium [element not in mmw table]
#   "NpF6":       (  7.19634,   1271.8225,   239.53,    54.85,    76.85),  # 270 neptunium hexafluoride [element not in mmw table]
    "O2":         (  6.83706,    339.2095,   268.70,  -218.80,  -118.57),  # 271 oxygen
    "O3":         (  7.15980,    618.6072,   255.87,  -193.00,   -12.15),  # 272 ozone
    "Os":         (  9.76546,  34751.0553,    35.61,  2600.00,  5400.00),  # 273 osmium
    "OsF5":       (  9.75828,   3432.7793,   273.23,    74.10,   225.90),  # 274 osmium fluoride
    "OsF6":       (  8.26547,   1525.6447,   235.83,   -89.00,    47.50),  # 275 osmium fluoride
    "OsOF5":      (  3.06277,     47.1535,     1.04,    31.65,    57.75),  # 276 osmium oxide pentafluoride
    "OsO4":       (  4.72185,    404.5748,    84.75,     3.20,   135.00),  # 277 osmium tetroxide
    "P":          (  7.00000,   2016.5981,   215.00,    69.00,   317.00),  # 278 phosphorus (white)
#   "P":          ( 11.39619,   5737.7187,   242.50,   181.85,   431.00),  # 279 phosphorus (red) [duplicate: red phosphorus; white phosphorus used instead]
    "PBrF2":      (  6.90420,    885.1200,   236.00,   -91.15,    10.85),  # 280 phosphorus bromide difluoride
    "PBr2F":      (  6.85800,   1210.3000,   226.00,   -18.15,   109.85),  # 281 phosphorus dibromide fluoride
    "PBr3":       (  7.35097,   1883.6717,   248.44,     7.80,   172.95),  # 282 phosphorus tribromide
    "PBr5":       ( 10.31441,   2746.6473,   259.11,   -19.00,   110.10),  # 283 phosphorus bromide
    "PClF2":      (  7.34733,   1041.3236,   280.39,  -119.50,    89.25),  # 284 phosphorus chloride difluoride
    "PCl2F":      (  7.11413,   1120.6064,   250.71,   -71.10,   189.85),  # 285 phosphorus dichloride fluoride
    "PCl2F3":     (  7.38018,   1175.5119,   257.82,   -79.80,   -22.80),  # 286 phosphorus dichloride trifluoride
    "PCl3":       (  6.87188,   1218.9462,   229.32,   -92.00,   101.00),  # 287 phosphorus chloride
    "PCl5":       ( 10.28827,   2973.0094,   241.50,   -12.00,   159.85),  # 288 phosphorus pentachloride
    "PF3":        (  7.06415,    651.8743,   257.33,  -152.00,    -1.95),  # 289 phosphorus fluoride
    "PF3O":       ( 10.28247,   1549.5473,   248.89,  -124.00,   -39.70),  # 290 phosphorus oxyfluoride
    "PF5":        (  9.61509,   1143.4832,   254.40,  -157.00,   -84.60),  # 291 phosphorus fluoride
    "PH3":        (  6.98375,    724.4708,   264.31,  -133.78,    51.60),  # 292 phosphine
    "PH4Br":      ( 10.09866,   2085.9953,   250.37,   -43.70,    38.30),  # 293 phosphonium bromide
    "PH4Cl":      (  9.08492,   1255.1903,   229.31,   -91.00,   -27.00),  # 294 phosphonium chloride
    "PH4I":       ( 10.23174,   2291.9882,   249.30,   -25.20,    62.50),  # 295 phosphonium iodide
    "POBrClF":    (  6.91410,   1214.0000,   222.00,   -10.15,   102.85),  # 296 phosphoryl bromide chloride fluoride
    "POBrCl2":    (  6.92410,   1411.0000,   213.00,    31.85,   162.85),  # 297 phosphoryl bromide dichloride
    "POBrF2":     (  7.10190,   1118.9000,   233.00,   -44.15,    51.85),  # 298 phosphoryl bromide difluoride
    "POBr3":      (  7.13796,   1697.4936,   207.04,    64.00,   191.70),  # 299 phosphorus oxybromide
    "POClF2":     (  6.92660,    946.9600,   231.00,   -76.15,    30.85),  # 300 phosphoryl chloride difluoride
    "POCl2F":     (  7.08465,   1201.8600,   233.00,   -40.15,    84.85),  # 301 phosphoryl dichloride fluoride
    "POCl3":      (  7.05109,   1390.0324,   227.82,     1.18,   105.50),  # 302 phosphorus oxychloride
    "PSBr3":      (  9.40218,   2659.4144,   232.93,    50.00,   175.00),  # 303 phosphorus thiobromide
    "PSCl3":      (  7.11462,   1501.0578,   229.54,   -36.20,   125.00),  # 304 phosphorus thiochloride
    "P2O3":       (  8.39335,   2597.3392,   298.17,    47.30,   173.00),  # 305 phosphorus trioxide
    "P2O5":       ( 10.74529,   6470.4817,   217.75,   230.00,   700.00),  # 306 phosphorus pentoxide
    "P2S5":       (  9.02309,   4959.9962,   293.52,   288.00,  1017.85),  # 307 phosphorus sulfide
    "P4O6":       (  8.35860,   2568.0727,   295.82,    39.70,   173.00),  # 308 phosphorous oxide
    "Pb":         (  7.46345,   9131.8353,   243.70,   434.85,  1749.00),  # 309 lead
    "PbBr2":      (  8.87105,   6987.2446,   274.44,   513.00,   892.00),  # 310 lead bromide
    "PbCl2":      (  7.64617,   5101.5890,   119.56,   547.00,   951.00),  # 311 lead chloride
    "PbF2":       (  7.62479,   6525.7272,    82.58,   861.00,  1293.00),  # 312 lead fluoride
    "PbI2":       (  6.76027,   3580.4550,    50.32,   479.00,   872.00),  # 313 lead iodide
    "PbO":        (  9.24183,  10820.2993,   227.03,   943.00,  1474.00),  # 314 lead oxide
#   "PbO":        (  9.33076,  11048.8015,   240.48,   724.00,  1471.00),  # 315 lead oxide (massicot) [duplicate: massicot; entry 314 used instead]
    "PbS":        (  8.84745,   7859.9264,    36.31,   852.00,  1281.00),  # 316 lead sulfide
    "Pd":         (  8.56766,  18403.4082,   273.13,  1300.00,  3400.00),  # 317 palladium
#   "Po":         (  7.18388,   5265.0104,   261.55,   500.00,  1050.00),  # 318 polonium [element not in mmw table]
#   "PoCl4":      (  9.33621,   4280.8340,   273.15,   300.60,   410.00),  # 319 polonium tetrachloride [element not in mmw table]
    "Pr":         (  6.78062,  14030.7383,    77.80,  1497.70,  3520.00),  # 320 praseodymium
    "PrCl3":      ( 10.44454,  15000.0000,   273.15,  1680.00,  1710.00),  # 321 praseodymium chloride
    "Pt":         (  8.58584,  22055.7478,    41.02,  2057.00,  3825.00),  # 322 platinum
#   "Pu":         (  8.08674,  18531.9774,   331.78,  1483.00,  3228.00),  # 323 plutonium [element not in mmw table]
#   "PuBr3":      (  8.15259,   8352.8990,    93.61,   681.00,  1476.85),  # 324 plutonium bromide [element not in mmw table]
#   "PuCl3":      (  7.22981,   8007.1303,    48.46,   760.00,  1776.85),  # 325 plutonium chloride [element not in mmw table]
#   "PuF3":       (  9.39122,  14145.3938,    49.57,  1396.00,  1726.85),  # 326 plutonium fluoride [element not in mmw table]
#   "Ra":         (  6.71596,   6622.4566,   200.78,   546.00,  1526.00),  # 327 radium [element not in mmw table]
    "Rb":         (  6.99992,   3926.1751,   265.16,   160.40,   688.00),  # 328 rubidium
    "RbBr":       (  8.13532,   8288.1158,   237.34,   781.00,  1340.00),  # 329 rubidium bromide
    "RbCl":       (  7.05736,   6116.0391,    74.38,   792.00,  1390.00),  # 330 rubidium chloride
    "RbF":        (  6.21681,   3433.9886,  -380.62,   921.00,  1410.00),  # 331 rubidium fluoride
    "RbI":        (  7.66815,   7039.8460,   170.51,   748.00,  1300.00),  # 332 rubidium iodide
    "Re":         (  8.97642,  34838.4703,   119.34,  3030.00,  5681.00),  # 333 rhenium
    "ReF5":       (  8.98978,   3010.3083,   271.47,    58.80,   221.30),  # 334 rhenium fluoride
    "ReF6":       (  8.47754,   1550.5300,   243.24,   -97.00,    33.80),  # 335 rhenium fluoride
    "ReF7":       ( 11.47517,   5114.2585,   521.55,    48.30,    71.85),  # 336 rhenium fluoride
    "ReOF4":      (  9.45214,   2533.4840,   213.84,     5.00,   171.70),  # 337 rhenium oxytetrafluoride
    "ReOF5":      (  7.90580,   1772.1930,   279.68,  -103.00,    73.00),  # 338 rhenium oxypentafluoride
    "ReO2F3":     ( 10.36000,   3437.0000,   273.15,    90.00,   181.85),  # 339 rhenium dioxytrifluoride
    "Re2O7":      ( 10.48807,   4082.7928,   176.70,   120.00,   500.00),  # 340 rhenium heptoxide
    "Rh":         (  9.68520,  26833.7062,   248.59,  2015.00,  3695.00),  # 341 rhodium
#   "Rn":         (  6.58011,    651.6732,   237.86,  -163.00,   103.85),  # 342 radon [element not in mmw table]
    "Ru":         (  9.33600,  26794.4200,    33.09,  2315.00,  4115.00),  # 343 ruthenium
    "RuF5":       (  7.56315,   2037.5032,   196.96,    49.60,   156.80),  # 344 ruthenium fluoride
    "S":          (  7.86630,   3464.9835,   250.34,   102.00,   444.00),  # 345 sulfur
    "SCl2":       (  8.36643,   1559.0062,   224.60,   -76.00,    59.60),  # 346 sulfur dichloride
    "SFN":        (  6.48180,    877.1000,   239.00,    -8.15,    29.85),  # 347 thiazyl fluoride
    "SF2N2":      (  6.78280,    901.0000,   242.00,   -83.15,     8.85),  # 348 dinitrogen sulfur difluoride
    "SF4":        (  7.02242,    882.8300,   253.61,  -100.00,   -40.45),  # 349 sulfur tetrafluoride
    "SF6":        (  6.87699,    676.2391,   233.12,  -168.00,    45.54),  # 350 sulfur hexafluoride
    "SOBr2":      (  8.17193,   2202.4959,   276.26,    -6.70,   140.00),  # 351 thionyl bromide
    "SOClF":      (  7.17310,   1100.1000,   244.00,   -65.15,    36.85),  # 352 thionyl chloride fluoride
    "SOCl2":      (  7.27960,   1415.8028,   246.26,  -101.15,    99.00),  # 353 thionyl chloride
    "SOF2":       (  6.28388,    569.2883,   211.09,   -98.73,   -43.80),  # 354 thionyl fluoride
    "SOF4":       (  7.07180,    840.3000,   249.00,  -113.15,   -25.15),  # 355 sulfur oxide tetrafluoride
    "SO2":        (  7.33311,   1013.4609,   237.65,   -75.48,   157.60),  # 356 sulfur dioxide
    "SO2BrF":     (  7.14280,   1155.0000,   231.00,   -45.15,    67.85),  # 357 sulfuryl bromide fluoride
    "SO2ClF":     (  6.52150,    793.7300,   210.70,   -72.15,    36.85),  # 358 sulfuryl chloride fluoride
    "SO2Cl2":     (  6.94482,   1194.8068,   224.60,   -51.15,   271.85),  # 359 sulfuryl chloride
    "SO2F2":      (  6.90700,    784.3000,   250.00,  -121.15,   -32.15),  # 360 sulfuryl fluoride
    "SO3":        (  6.85638,    754.8178,   145.11,    16.80,   217.70),  # 361 sulfur trioxide
    "S2Br2":      (  7.17834,   1664.6387,   185.93,    -7.00,   200.90),  # 362 sulfur bromide
    "S2Cl2":      (  7.04178,   1466.8804,   215.53,    -7.40,   137.00),  # 363 sulfur monochloride
    "S2F2":       (  6.68410,    628.0000,   256.00,  -120.15,   -77.15),  # 364 disulfur difluoride
    "S2F10":      (  7.06760,   1100.6000,   234.00,   -47.15,    48.85),  # 365 sulfur decafluoride
    "S2O2F10":    (  6.87410,   1110.0000,   229.00,   -34.15,    70.85),  # 366 thiosulfuryl decafluoride
    "S2O5ClF":    (  7.01560,   1257.4000,   204.00,     4.85,   130.85),  # 367 pyrosulfuryl chloride fluoride
    "S2O5Cl2":    (  7.01910,   1460.0000,   202.00,    47.85,   177.85),  # 368 pyrosulfuryl chloride
    "S2O5F2":     (  6.88110,   1120.0000,   229.00,   -33.15,    72.85),  # 369 pyrosulfuryl difluoride
    "S2O5F4":     (  6.88510,   1140.0000,   227.00,   -27.15,    79.85),  # 370 disulfur pentoxide tetrafluoride
    "S3O8F2":     (  6.99010,   1360.0000,   211.00,    22.85,   145.85),  # 371 trisulfur octoxide difluoride
    "Sb":         (  5.88125,   4984.6752,    74.32,   534.00,  1587.00),  # 372 antimony
    "SbBr3":      (  8.00960,   2654.1840,   237.51,    93.90,   280.00),  # 373 antimony tribromide
    "SbCl3":      (  6.91949,   1709.1498,   202.95,    73.40,   520.85),  # 374 antimony trichloride
    "SbCl5":      (  7.76515,   2046.6707,   240.86,    22.70,   114.10),  # 375 antimony pentachloride
    "SbH3":       (  7.20112,   1062.6666,   263.97,   -95.28,   167.20),  # 376 stibine
    "SbI3":       (  7.55355,   2915.0705,   222.85,   163.60,   401.00),  # 377 antimony triiodide
    "Sb2O3":      (  4.06652,   1464.9462,  -189.49,   400.00,  1550.00),  # 378 antimony oxide
    "Sc":         (  7.72648,  13961.2550,    45.18,  1300.00,  3200.00),  # 379 scandium
    "ScBr3":      ( 14.35000,  13780.0000,   273.15,   768.85,   926.85),  # 380 scandium bromide
    "ScCl3":      ( 14.37000,  14200.0000,   273.15,   791.85,   959.85),  # 381 scandium chloride
    "ScI3":       ( 14.17000,  13340.0000,   273.15,   736.85,   906.85),  # 382 scandium iodide
    "Se":         (  7.05403,   3441.9628,   139.77,   227.00,   685.00),  # 383 selenium (gray)
    "SeCS":       (  6.89093,   1274.4371,   232.21,   -47.30,    85.60),  # 384 carbonyl selenide sulfide
    "SeCl4":      ( 10.79338,   3477.8728,   248.14,    74.00,   191.40),  # 385 selenium tetrachloride
    "SeF4":       (  7.61816,   1471.6008,   204.64,    13.60,   106.00),  # 386 selenium tetrafluoride
    "SeF6":       (  6.99021,    727.3425,   222.79,  -153.00,    72.35),  # 387 selenium hexafluoride
    "SeOCl2":     (  7.90109,   1958.4938,   213.12,    34.80,   177.00),  # 388 selenium oxychloride
    "SeOF2":      (  7.42010,   1380.0000,   178.00,    36.85,   152.85),  # 389 selenium oxyfluoride
    "SeO2":       ( 18.01392,  11598.1403,   451.41,   115.00,   335.00),  # 390 selenium dioxide
    "Si":         (  8.71307,  20614.1581,   267.17,  1634.85,  3264.85),  # 391 silicon
    "SiBrCl2F":   (  7.51839,   1472.7873,   282.18,   -86.50,   201.15),  # 392 silane bromodichlorofluoro
    "SiBr2ClF":   (  7.26677,   1379.5487,   255.04,   -65.20,   242.77),  # 393 silane dibromochlorofluoro
    "SiCCl3N":    (  7.75100,   1687.0000,   273.15,   -46.15,    73.24),  # 394 silane trichloro-carbonitrile
    "SiCCl3NS":   (  7.99850,   2060.0000,   273.15,    66.85,   129.38),  # 395 silane trichloroisothio-cyanato
    "SiCH2Cl4O":  (  8.39500,   2159.0000,   273.15,    -0.15,   118.39),  # 396 silane trichloro(chloro-methoxy)
    "SiCH3Cl2F":  (  7.74530,   1472.3000,   273.15,   -30.15,   -27.15),  # 397 silane dichlorofluoro-methyl
    "SiCH3Cl3":   (  7.18818,   1370.6632,   251.81,   -77.80,   243.85),  # 398 silane methyltrichloro
    "SiCH3F3":    (  7.97520,   1238.1000,   273.15,   -60.15,   -57.15),  # 399 silane trifluoromethyl
    "SiCH4ClF":   (  7.60660,   1296.6000,   273.15,   -75.15,     1.22),  # 400 silane chlorofluoro-methyl
    "SiCH4Cl2":   (  7.23413,   1272.2113,   250.69,   -90.60,   209.85),  # 401 silane methyl dichloro
    "SiCH4F2":    (  7.76550,   1160.4000,   273.15,   -75.15,   -35.59),  # 402 silane difluoromethyl
    "SiCH5Cl":    (  7.04725,   1054.7751,   244.46,  -134.10,   168.85),  # 403 silane methyl chloro
#   "SiCH5Cl":    (  8.09302,   1719.2800,   298.60,   -27.15,    31.26),  # 404 silane (chloromethyl) [duplicate: isomer; entry 403 used instead]
    "SiCH6":      (  7.10688,    873.3901,   263.57,  -156.81,    79.35),  # 405 silane methyl
    "SiCH6O":     (  3.98012,    185.6317,   149.98,   -90.20,    18.88),  # 406 methyl silyl ether
    "SiC2H8":     (  6.73312,    853.3608,   241.12,  -150.22,   128.85),  # 407 silane dimethyl
    "SiC3H10":    (  6.90193,   1030.4100,   249.55,  -135.89,   158.85),  # 408 silane trimethyl
    "SiC4H12":    (  6.87190,   1060.6645,   239.11,   -99.08,   177.25),  # 409 silane tetramethyl
    "SiClF3":     (  7.06315,    761.3693,   252.04,  -144.00,    35.68),  # 410 silane chlorotrifluoro
    "SiCl2F2":    (  7.33305,   1054.7676,   268.91,  -124.70,    94.20),  # 411 silane dichlorodifluoro
    "SiCl3F":     (  7.15967,   1116.6266,   248.71,   -92.60,   161.70),  # 412 silane trichlorofluoro
    "SiCl4":      (  7.00381,   1205.9453,   235.64,   -68.85,   233.85),  # 413 silicon tetrachloride
    "SiF4":       (  6.72479,    459.3741,   215.21,  -156.00,   -14.15),  # 414 silicon tetrafluoride
    "SiHBr3":     (  6.97115,   1382.4891,   228.99,   -30.50,   344.35),  # 415 silane tribromo
    "SiHCl3":     (  7.10586,   1191.2883,   250.11,  -128.20,   205.85),  # 416 silane trichloro
    "SiHF3":      (  6.50485,    466.9478,   223.85,  -152.00,     3.50),  # 417 silane trifluoro
    "SiH2Br2":    (  7.05269,   1295.7656,   244.60,   -60.90,   286.09),  # 418 silane dibromo
    "SiH2Cl2":    (  6.98538,   1032.4344,   243.23,  -122.00,   175.85),  # 419 silane dichloro
    "SiH2F2":     (  6.68670,    609.3862,   237.92,  -146.70,    37.60),  # 420 silane difluoro
    "SiH2I2":     (  7.36622,   2014.1634,   299.05,     3.80,   410.66),  # 421 silane diiodo
    "SiH3Br":     (  7.24230,   1152.7870,   261.91,   -85.70,   182.00),  # 422 silane monobromo
    "SiH3Cl":     (  6.76710,    797.8130,   235.69,  -117.80,   123.64),  # 423 silane monochloro
    "SiH3F":      (  6.29474,    412.9780,   218.97,  -153.00,    12.72),  # 424 silane monofluoro
    "SiH3I":      (  7.04589,   1204.6095,   243.82,   -53.00,   251.44),  # 425 silane iodo
    "SiH4":       (  7.12525,    704.8874,   278.22,  -185.00,    -3.45),  # 426 silane
    "SiI4":       (  6.97269,   1909.4768,   201.81,   120.50,   256.85),  # 427 silane tetraiodo
    "SiO2":       (  9.87451,  16418.2297,   117.52,  1609.85,  2230.05),  # 428 silicon dioxide (alpha quartz)
    "Si2CH3Cl5":  (  7.91847,   2186.1100,   281.29,    34.85,    37.85),  # 429 silane pentachloromethyl-di
    "Si2CH7Cl":   (  7.83841,   1674.1300,   273.15,    -0.15,    64.54),  # 430 silane [(chlorosilyl)methyl]-
#   "Si2CH7Cl":   (  8.09282,   1875.0700,   273.15,   -23.15,    86.61),  # 431 silane (chloromethyl)-di [duplicate: isomer; entry 430 used instead]
    "Si2C2H4Cl6": (  7.82031,   2242.0000,   251.00,    90.85,   202.89),  # 432 1,2-bis(trichlorosilyl)ethane
#   "Si2C2H4Cl6": (  8.10890,   2425.0000,   273.15,    54.85,   190.69),  # 433 1,1,1,3,3,3-hexachloro-2-methyl-1,3-disilapropane [duplicate: isomer; entry 432 used instead]
    "Si2C2H6Cl4": (  7.14487,   1680.7900,   236.84,    26.85,   157.34),  # 434 silane 1,2,2,2-tetrachloro-1,1-dimethyldi
    "Si2C2H9FO":  (  6.32253,   1162.7400,   274.47,   -78.15,    63.37),  # 435 1-fluoro-1,3-dimethyldisiloxane
    "Si2Cl3F3O":  (  9.03398,   2311.2185,   332.87,   -38.00,    42.74),  # 436 disiloxane, trichlorotrifluoro
    "Si2Cl6":     (  7.36013,   1615.6257,   215.69,     4.00,   145.00),  # 437 disilane hexachloro
    "Si2Cl6O":    (  7.12363,   1465.7767,   210.43,    -5.00,   135.60),  # 438 disiloxane, hexachloro
    "Si2F6":      ( 12.22960,   2464.6767,   282.54,   -81.00,   -18.90),  # 439 disilane hexafluoro
    "Si2H6":      (  7.00941,   1015.4056,   260.09,  -129.30,   158.85),  # 440 disilane
    "Si2H6O":     (  7.42748,   1137.4567,   265.37,  -112.50,   -15.20),  # 441 disiloxane
    "Si3Cl8":     (  8.65954,   2872.1180,   285.62,    46.30,   211.40),  # 442 octachlorotrisilane
    "Si3H8":      (  7.22385,   1329.1653,   253.14,   -68.90,    52.90),  # 443 trisilane
    "Si3H9N":     (  7.15693,   1247.8657,   243.12,   -68.70,    48.70),  # 444 trisilazane
    "Si4H10":     (  6.54087,   1122.5646,   199.31,   -27.70,   107.40),  # 445 tetrasilane
    "Sm":         (  6.72424,   7242.8476,    90.48,   700.00,  2000.00),  # 446 samarium
    "Sn":         (  8.54915,  16655.8123,   336.40,  1150.00,  2800.00),  # 447 tin (white)
    "SnBr4":      (  7.14579,   1846.6483,   227.98,    58.30,   205.00),  # 448 stannic bromide
    "SnCl2":      (  7.51627,   3705.3738,   176.36,   316.00,   623.00),  # 449 stannous chloride
    "SnCl4":      (  6.85720,   1294.1344,   211.30,   -22.70,   114.15),  # 450 stannic chloride
    "SnH4":       (  7.15934,    937.8599,   271.00,  -140.00,   -51.80),  # 451 stannane
    "SnI4":       (  6.06085,   1628.2894,   147.69,   156.00,   364.35),  # 452 stannic iodide
    "SnS":        (  9.97000,  10460.0000,   273.15,   726.85,  1026.85),  # 453 tin sulfide
    "Sr":         (  6.83233,   6128.6851,   168.97,   523.00,  1382.00),  # 454 strontium
    "SrF2":       (  8.06789,  12904.6141,     3.46,  1477.00,  2426.85),  # 455 strontium fluoride
    "SrO":        (  7.51441,   9499.3367,  -803.52,  1750.00,  2700.00),  # 456 strontium oxide
    "Ta":         (  8.73303,  31263.6507,  -115.81,  3024.00,  5458.00),  # 457 tantalum
    "TaBr5":      (  8.14000,   3260.0000,   273.15,   265.00,   336.85),  # 458 tantalum bromide
    "TaCl5":      ( 11.98340,   4654.0000,   273.15,   188.13,   238.13),  # 459 tantalum chloride
    "TaF5":       (  6.45981,   1797.9105,   273.15,   119.00,   249.20),  # 460 tantalum fluoride
    "Tb":         (  7.41799,  14824.6184,    37.37,  1516.10,  3230.00),  # 461 terbium
#   "Tc":         ( 26.08539, 236826.8849,  5941.04,  2400.00,  4400.00),  # 462 technetium [element not in mmw table]
#   "Tc2O7":      (  9.00000,   3570.0000,   273.15,   118.85,   306.85),  # 463 technetium heptoxide [element not in mmw table]
    "Te":         (  6.63850,   4084.3774,    98.94,   502.00,   988.00),  # 464 tellurium
    "TeCl2":      (  8.51000,   3350.0000,   273.15,   208.00,   316.85),  # 465 tellurium dichloride
    "TeCl4":      (  8.59582,   3538.9650,   232.24,   233.00,   387.00),  # 466 tellurium tetrachloride
    "TeF4":       (  5.64000,   1787.0000,   273.15,   129.00,   196.85),  # 467 tellurium tetrafluoride
    "TeF6":       (  8.65918,   1258.5289,   256.70,  -111.30,   -38.90),  # 468 tellurium hexafluoride
    "TeO2":       ( 12.76126,  15000.0000,   273.15,  1215.00,  1245.00),  # 469 tellurium dioxide
    "Te2F10":     (  6.90180,   1150.0000,   227.00,   -33.15,    89.85),  # 470 tellurium decafluoride
    "Th":         (  8.08673,  25785.3686,   165.09,  2360.00,  4788.00),  # 471 thorium
    "ThBr4":      (  9.56000,   7550.0000,   273.15,   679.00,   852.85),  # 472 thorium bromide
    "ThCl4":      (  9.57000,   7980.0000,   273.15,   770.00,   912.85),  # 473 thorium chloride
    "ThI4":       (  9.09000,   6890.0000,   273.15,   582.85,   833.85),  # 474 thorium iodide
    "ThO2":       (  9.64285,  31600.0000,   273.15,  4350.00,  4400.00),  # 475 thorium oxide
    "Ti":         (  8.90223,  20948.9918,   190.76,  1709.00,  3285.00),  # 476 titanium
    "TiBr4":      (  6.97018,   1791.8855,   205.54,    39.00,   228.85),  # 477 titanium bromide
    "TiCl2":      (  8.35740,  10003.3561,   174.30,  1035.00,  1426.85),  # 478 titanium chloride
    "TiCl4":      (  7.29500,   1668.8244,   242.21,   -24.10,   364.85),  # 479 titanium tetrachloride
    "TiF4":       ( 10.87761,   4455.4168,   273.15,   234.00,   284.00),  # 480 titanium fluoride
    "TiI4":       (  7.57600,   3054.0000,   273.15,   156.85,   369.85),  # 481 titanium iodide
    "Tl":         (  7.58000,   8060.0496,   229.50,   609.00,  1485.00),  # 482 thallium
    "TlBr":       (  7.32435,   4432.7693,   178.58,   490.00,   819.00),  # 483 thallous bromide
    "TlCl":       ( -5.69389,   6065.8547, -1427.41,   504.00,   720.00),  # 484 thallium chloride
    "TlF":        (  9.88714,   6902.2858,   273.15,   662.00,   712.00),  # 485 thallium fluoride
    "TlI":        (  7.34352,   4374.7071,   156.28,   440.00,   824.00),  # 486 thallium iodide
    "Tm":         (  7.82243,  10435.3130,   205.05,   844.00,  1297.00),  # 487 thulium
    "U":          (  7.77000,  23272.8375,   600.00,  2052.00,  4129.00),  # 488 uranium
    "UBr3":       (  9.06095,  10506.7518,   178.05,   727.00,   926.85),  # 489 uranium bromide
    "UBr4":       (  8.12161,   4834.9448,   143.94,   519.00,   766.85),  # 490 uranium bromide
    "UCl4":       (  9.01541,   5767.4399,   147.57,   590.00,   770.00),  # 491 uranium chloride
    "UCl6":       (  7.87988,   4000.0000,   273.15,   477.00,   527.00),  # 492 uranium chloride
    "UF4":        (  8.93038,   9333.8078,    95.95,  1036.00,  1400.00),  # 493 uranium fluoride
    "UF6":        (  8.47645,   1568.8742,   223.87,   -38.80,    56.50),  # 494 uranium fluoride
    "UI4":        (  8.81552,   5446.7685,   159.45,   506.00,   726.85),  # 495 uranium iodide
    "V":          (  9.18000,  22303.0522,   133.69,  1828.00,  3406.85),  # 496 vanadium
    "VCl4":       (  7.32955,   1823.3734,   258.01,   -25.70,   423.85),  # 497 vanadium tetrachloride
    "VF5":        ( 10.43000,   2423.0000,   273.15,    19.50,    45.00),  # 498 vanadium fluoride
    "VOCl3":      (  7.02483,   1518.9432,   239.69,   -79.50,   126.85),  # 499 vanadium oxytrichloride
    "V2O5":       (  5.05000,   7100.0000,   273.15,   670.00,  1226.85),  # 500 vanadium oxide
    "W":          (  9.98307,  40387.9650,   131.64,  3204.00,  5555.00),  # 501 tungsten
    "WCl4":       (  8.19500,   3253.0000,   273.15,   281.85,   326.85),  # 502 tungsten chloride
    "WCl5":       (  7.72000,   2760.0000,   273.15,   242.00,   280.00),  # 503 tungsten chloride
    "WCl6":       (  7.87000,   3050.0000,   273.15,   275.00,   340.00),  # 504 tungsten chloride
    "WF6":        (  7.07104,    907.7240,   199.53,   -71.60,    17.10),  # 505 tungsten fluoride
    "WOF4":       (  9.32296,   2709.3726,   234.67,     2.00,   185.90),  # 506 tungsten oxytetrafluoride
    "WO3":        ( 14.53875,  24600.0000,   273.15,  1787.00,  1837.00),  # 507 tungsten oxide
    "Xe":         (  7.02795,    701.1036,   277.06,  -111.75,    16.55),  # 508 xenon
    "XeF2":       (  6.96427,   1073.3758,   148.51,     2.90,   357.85),  # 509 xenon difluoride
    "XeF4":       ( 10.91390,   3095.0600,   269.56,    35.85,   117.85),  # 510 xenon tetrafluoride
    "Y":          (  8.00530,  17991.8060,   165.95,  1610.10,  3345.00),  # 511 yttrium
    "Yb":         (  7.28084,   6060.1722,   181.30,   463.00,  1196.00),  # 512 ytterbium
    "Zn":         (  7.57409,   5241.2423,   203.39,   337.00,   912.00),  # 513 zinc
    "ZnBr2":      (  8.49304,   4782.6046,   185.14,   397.85,   561.85),  # 514 zinc bromide
    "ZnCl2":      (  8.50869,   5054.5539,   166.13,   428.00,   732.00),  # 515 zinc chloride
    "ZnF2":       (  6.56144,   4368.9367,  -312.99,   970.00,  1500.00),  # 516 zinc fluoride
    "ZnI2":       ( 18.01392,  19725.9709,   678.50,   301.00,   625.00),  # 517 zinc iodide
    "ZnO":        (  9.97977,  15781.3538,   273.15,  1499.95,  1949.95),  # 518 zinc oxide
    "ZnS":        ( 12.46830,  13980.0000,   273.15,  1135.00,  1185.00),  # 519 zinc sulfide (wurtzite)
    "ZnSO4":      (  8.74032,   2215.7753,   273.15,    20.00,   105.00),  # 520 zinc sulfate
    "Zr":         (  9.12202,  28773.5320,   201.25,  2366.00,  4409.00),  # 521 zirconium
    "ZrBr4":      ( 10.28269,   4025.9186,   183.91,   207.00,   360.00),  # 522 zirconium bromide
    "ZrCl4":      ( 11.44848,   4801.8756,   229.46,   190.00,   331.00),  # 523 zirconium chloride
    "ZrF4":       ( 15.61264,  15000.0000,   273.15,   875.00,   905.00),  # 524 zirconium fluoride
    "ZrI4":       ( 11.20709,   5397.4287,   217.24,   264.00,   431.00),  # 525 zirconium iodide
}
yaws_cite = "C.L. Yaws and M.A. Satyro, The Yaws Handbook of Vapor Pressure (Second Edition), Chapter 2, Elsevier 2015, doi:10.1016/B978-0-12-802999-2.00002-7"

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

# Piecewise fit of the form ln(P/bar) = a - b/T
def lnfit_psat(t:float, gas:str):
    fit = gases[gas]["fit_ln"]

    A = []
    for f in fit:
        t_min = f[0]
        if t > t_min:
            A = [f[1],f[2]]

    if len(A) == 0:
        raise Exception("Temperature out of range")

    return np.exp( A[0] - A[1]/t ) * bar  # Pa


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

def yaws_load_all():
    # add Yaws & Satyro (2015) fits for species not already covered above or by Woitke+2017
    added = []
    for gas, (A, B, C, tmin, tmax) in yaws_inorganic.items():
        if (gas in gases) or (gas in woitke_gases) or (gas == "H2O"):
            continue
        gases[gas] = {
            "cite": [yaws_cite],
            "fit":  [   # convert to P in bar and T in kelvin
                        [tmin + 273.15, A - np.log10(bar/mmHg), B, C - 273.15]
                    ],
            # Yaws range of validity, not necessarily the true critical and triple points
            "T_crit": tmax + 273.15,
            "T_trip": tmin + 273.15,
            "prefer": "antoine",
        }
        added.append(gas)
    return added

yaws_gases = yaws_load_all()
print("Loaded Yaws & Satyro (2015) Antoine coefficients for %d gases"%len(yaws_gases))

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
        