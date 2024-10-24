# Shmoscar: Core module of SHMOSCAR
# This is the place I chuck my global variables and most heavily-used functions

import re # impt: read string of numbers
import math
import numpy as np
import pandas as pd
import itertools
from operator import add
from copy import deepcopy

# Global Variables
atom_name_index = 5
atom_number_index = 6

# Text elements
longspace = "    "
ls = longspace # just to make code sexier
lls = ls + "   " # extra long longspace
horizont = "=============================================\n"
flpr = "{:11f}" # float precision string (adjust here to taste)

# For people to check if input atomic symbol is a real element
periodic_table = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']

def readfile(in_filename, verbose = True):
    # Read input file
    infile = open(in_filename, "r")
    data = infile.readlines()
    infile.close()
    if verbose: print("Reading from file: " + in_filename)
    return data

def writefile(data, out_filename, verbose = True):
    # Write file data (data) to output file (out_filename)
    outfile = open(out_filename, "w")
    outfile.writelines(data)
    outfile.close()
    if verbose: print("Writing to file: " + out_filename)
    
def rename(data, newname):
    # Replace description of POSCAR file
    # Replace header
    data[0] = newname + "\n"
    
def basis(data):
    # Read lattice vectors
    m1 = re.findall(r"-?\d+\.\d+", data[1].strip())
    mu = float(m1[0])
    B = np.zeros((3,3))
    for i in range(3):
        basis = np.array(re.findall(r"-?\d+\.\d+", data[i+2].strip()))
        B[i] = mu * basis.astype(float)
    return B
    
def printvaspdata(data):
    print(">>>>>>>>>>>> START VASPFILE >>>>>>>>>>>>\n")
    print("".join(data))
    print(">>>>>>>>>>>>> END VASPFILE >>>>>>>>>>>>>\n")