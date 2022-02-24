# Cartesian: Convert Direct to Cartesian coordinates, and vice versa 

# Andy Paul Chen, Tuesday 21 April 2020, Cleveland Heights, Ohio

# 18 February 2022: updated to account for line 1 multiplier of lattice vectors

from Seldyn import * # selective dynamics package
ls = longspace

def data_isCart(data):
    dcindex = 8 if data_isSeldyn(data) else 7 # Index of Direct/Cartesian line
    if data[dcindex][0].upper() == 'C' or data[dcindex][0].upper() == 'K':
        return True
    
def isCart(in_filename):
    # Read input file
    data = readfile(in_filename)
    if data_isCart(data):
        print("Coordinates form: Cartesian")
    else:
        print("Coordinates form: Fractional")
        
def data_switchCart(data):
    # Switches coordinates between Cartesian and Direct!
    dcindex = 8 if data_isSeldyn(data) else 7 # Index of Direct/Cartesian line
    num = re.findall(r'\d+', data[atom_number_index].strip())
    snum = sum(list(map(int,num)))
    data = data[:(snum + dcindex + 1)] # whittle the tail
    
    # Read lattice vectors
    m1 = re.findall(r"-?\d+\.\d+", data[1].strip())
    mu = float(m1[0])
    a = re.findall(r"-?\d+\.\d+", data[2].strip())
    b = re.findall(r"-?\d+\.\d+", data[3].strip())
    c = re.findall(r"-?\d+\.\d+", data[4].strip())
    
    # Conversion operation
    if data_isCart(data):
        print("Converting: Cartesian > Fractional (not supported in this version)")
    else:
        print("Converting: Fractional > Cartesian")
        for line in range(len(data)):
            if line > dcindex:
                dirs = re.findall(r"-?\d+\.\d+", data[line].strip())
                if data_isSeldyn(data):
                    flags = re.findall(r"\w", data[line].strip())[-3:]
                carts = [0,0,0]
                for i in range(3):
                    carts[i] += mu*float(a[i])*float(dirs[0]) + mu*float(b[i])*float(dirs[1]) + mu*float(c[i])*float(dirs[2])
                data[line] = ls + "{:11f}".format(carts[0]) + ls + "{:11f}".format(carts[1]) + ls + "{:11f}".format(carts[2])
                if data_isSeldyn(data):
                    data[line] += longspace + flags[0] + " " + flags[1] + " " + flags[2]
                data[line] += "\n"
    
    # Switch the Direct/Cartesian line
    lineindex = 8 if data_isSeldyn(data) else 7
    if data_isCart(data):
        data[lineindex] = "Direct\n"
    else:
        data[lineindex] = "Cartesian\n"
    return data
    
def switchCart(in_filename, out_filename):
    # Read input file
    data = readfile(in_filename)
    
    data = data_switchCart(data)
    
    # Write output file
    writefile(data, out_filename)