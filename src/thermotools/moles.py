import numpy as np
import os
from thermotools import get_inpdir, empty_dir

def formula_from_path(f:str):
    return f.split("/")[-1].split(".")[0]

def read_elements():
    # Read table of molecular weights [kg/mol]
    fpath = os.path.join(get_inpdir(),"mmw","web","elements.txt")

    with open(fpath,'r') as hdl:
        lines = hdl.readlines()

    elem = {}
    for l in lines[2:]:
        s = l.split()

        # key 
        k = s[1].strip()

        # value 
        v = s[3].strip()
        v = v.split("(")[0]
        if "[" in v:
            continue 
        v = float(v)

        # store 
        elem[k] = v

    head = "Symbol , Mass [kg mol-1]"
    X = []
    X_dict = {}
    for k in elem.keys():
        mu = elem[k]*1.0e-3
        k = str(k)
        X.append([k,"%.9e"%mu])
        X_dict[k] = mu

    datdir = os.path.join(get_inpdir(),"mmw","dat")
    empty_dir(datdir)
    np.savetxt(os.path.join(datdir,"elements.csv"),X,fmt=["%s","%s"],delimiter=',',header=head)

    return X_dict


# Count atoms in a molecule 
def count_atoms(m:str):
        # Setup 
        out = {}
        nchar = len(m)
        i = 0
        e = ""
        count = -1
        last = False

        # Loop through string
        while i < nchar:
            last = (i == nchar-1)

            # new element 
            # print(m[i])
            if m[i].isupper():
                count = 0
                e = str(m[i])
                if (not last) and (m[i+1].islower()):  # two letter element name 
                    e = e+str(m[i+1])
                    i += 1

            last = (i == nchar-1)

            # get count 
            if count == 0:   # expecting number
                # number of atoms 
                if last or m[i+1].isalpha(): # got letter => count=1
                    count = 1
                else:
                    count = int(m[i+1]) 

                # repeated element 
                if e in out.keys():
                    out[e] += count 
                else:
                    out[e] = count 

                # reset 
                e = ""
                count = -1 
            i += 1
        
        return out 

# Calculate mmw from formula
def mmw_from_formula(m:str, elem_table:dict):

    # get atoms 
    atoms = count_atoms(m)

    # add up atoms 
    mmw = 0.0
    for k in atoms.keys():
        mmw += elem_table[k]*atoms[k]

    return mmw

