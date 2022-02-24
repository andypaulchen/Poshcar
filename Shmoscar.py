# Shmoscar: Core module of SHMOSCAR
# This is the place I chuck my global variables and most heavily-used functions

import re # impt: read string of numbers
import math
import numpy as np

# Global Variables
atom_name_index = 5
atom_number_index = 6
longspace = "    "
horizont = "============================================="

# For people to check if input atomic symbol is a real element
Periodic_Table = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']

def readfile(in_filename):
    # Read input file
    infile = open(in_filename, "r")
    data = infile.readlines()
    infile.close()
    print("Reading from file: " + in_filename)
    return data

def writefile(data, out_filename):
    # Write file data (data) to output file (out_filename)
    outfile = open(out_filename, "w")
    outfile.writelines(data)
    outfile.close()
    print("Writing to file: " + out_filename)
    
def rename(in_filename, newname):
    # Replace description of 
    # Read input file
    data = readfile(in_filename)
    
    # Replace header
    data[0] = newname + "\n"
    
    # Write output file
    writefile(data, in_filename)