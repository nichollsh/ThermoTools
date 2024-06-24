
import numpy as np 
import matplotlib.pyplot as plt 
import os, glob, time
import netCDF4 as nc

today = str(time.strftime("%Y-%m-%d"))

import moles,  phase,   water


tmin = 0.5
dt   = 0.5

print("START")

print("Clear old files")
for fold in ["sat","lv","compiled"]:
    for f in glob.glob(fold+"/dat/*"):
        os.remove(f)

print("Read elements")
elem_table = moles.read_elements()



# --------------------------------------------------
# LATENT HEAT 
# --------------------------------------------------

def write_csv_lat(gas,ts,hs):
    X = np.array([ts,hs]).T 
    head = "T [K] , delta-h [J kg-1]"
    fpath = "lv/dat/%s.csv"%gas 

    np.savetxt(fpath, X, fmt="%.9e", delimiter=',', header=head)


print("\nLatent heat curve")
# Coker (2007) gases (vapourisation only)
data_vap = np.loadtxt("lv/web/coker2007_appendixC.csv", delimiter=',', dtype=str, skiprows=1)
coker_gases = []
for row in data_vap:

    # get values
    formula = row[0]
    name    = row[1]
    fit_A   = float(row[2])
    fit_T   = float(row[3])
    fit_n   = float(row[4])
    Tmin    = float(row[5])
    Tmax    = float(row[6])

    coker_gases.append(formula)
    print("\t",formula, "coker2007")

    # define fit 
    def _fit(t:float):
        t = max(t, Tmin)
        t = min(t, Tmax)
        out = fit_A * (1 - t/fit_T)**fit_n # kJ/mol
        out = out * 1e3 / moles.mmw_from_formula(formula, elem_table)  # J/kg
        return out 
    
    t_arr = np.arange(tmin, Tmax+dt, dt)
    h_arr = [_fit(t) for t in t_arr]
    write_csv_lat(formula, t_arr, h_arr)


# Gases without h_vap data
no_vap = []
for k in phase.gases.keys():
    gas = phase.gases[k]

    if k in coker_gases:
        continue 
    if k == "H2O":
        continue 

    no_vap.append(k)
    print("\t",k, "no_data")
    write_csv_lat(k, [0.0, gas["T_crit"]], [0.0,0.0])

# Water
print("\t H2O W&P")
t_arr = np.arange(tmin, water.T_crit, dt)[:-1]
t_arr = np.concatenate((t_arr, [water.T_crit]))
h_arr = [water.delta_both(t) for t in t_arr]
write_csv_lat("H2O",t_arr, h_arr)

# Make plot 
print("\t Plot")
fig,ax = plt.subplots(1,1)
files = glob.glob("lv/dat/*.csv")
for i,f in enumerate(files):
    form = moles.formula_from_path(f)
    data = np.loadtxt(f, delimiter=',').T

    ls='solid'
    if i > 9:
        ls='dashed'
    if i > 18:
        ls='dotted'

    ax.plot(data[0], data[1], label=form, ls=ls)

ax.set(xlabel="Temperature [K]", ylabel="Enthalpy of phase change [J/kg]")
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncols=2)
fig.savefig("lv/all.png", dpi=220, bbox_inches='tight')
plt.close("all")


# --------------------------------------------------
# SATURATION CURVE
# --------------------------------------------------

print("\nSaturation curve...")

def write_csv_sat(gas,ts,ps):
    X = np.array([ts,ps]).T 
    head = "T [K] , P_sat [Pa]"
    fpath = "sat/dat/%s_sat.csv"%gas 

    np.savetxt(fpath, X, fmt="%.9e", delimiter=',', header=head)

def write_tripcrit(gas, t_trip, t_crit):

    X = [t_trip]
    head = "T_trip [K]"
    fpath = "sat/dat/%s_trip.csv"%gas 
    np.savetxt(fpath, X, fmt="%.9e", delimiter=',', header=head)

    X = [t_crit]
    head = "T_crit [K]"
    fpath = "sat/dat/%s_crit.csv"%gas 
    np.savetxt(fpath, X, fmt="%.9e", delimiter=',', header=head)

# Main gases
for gas in coker_gases:

    if gas == "H2O":
        continue 
    if gas not in phase.gases.keys():
        print("\t",gas,"SKIPPED")
        continue 

    print("\t",gas)
    gdict = phase.gases[gas]
    mmw =  moles.mmw_from_formula(gas, elem_table)

    t_arr = np.arange(tmin, gdict["T_crit"], dt)
    t_arr = np.concatenate((t_arr, [gdict["T_crit"]]))
    p_arr = [phase.cc_psat(t, gas, mmw) for t in t_arr]
    write_csv_sat(gas, t_arr, p_arr)
    write_tripcrit(gas, gdict["T_trip"], gdict["T_crit"])

# Water
print("\t","H2O")
t_arr = np.arange(tmin, water.T_crit, dt)
t_arr = np.concatenate((t_arr, [water.T_crit]))
arr_p = [water.psat_both(t) for t in t_arr]
write_csv_sat("H2O",t_arr, arr_p)
write_tripcrit("H2O", water.T_trip, water.T_crit)

# Plot 
print("\t","Plot")
fig,ax = plt.subplots(1,1)
files = glob.glob("sat/dat/*_sat.csv")
for i,f in enumerate(files):
    form = moles.formula_from_path(f).split("_")[0]
    data = np.loadtxt(f, delimiter=',').T

    ls='solid'
    if i > 9:
        ls='dashed'
    if i > 18:
        ls='dotted'

    l = ax.plot(data[0], data[1]*1e-5, label=form, ls=ls, zorder=3)[0]

ax.set(xlabel="Temperature [K]", ylabel="Saturation pressure [bar]")
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncols=2)
ax.set_yscale("log")
ax.set_ylim(bottom=1e-10, top=1e4)
fig.savefig("sat/all.png", dpi=220, bbox_inches='tight')
plt.close("all")




# --------------------------------------------------
# COMPILE TO NETCDF
# --------------------------------------------------

files = glob.glob("lv/dat/*.csv")
gases = [moles.formula_from_path(f) for f in files]
elem_table = moles.read_elements()
print(gases)