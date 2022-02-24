# VacuumQD: Create a QD surrounded by vacuum from a single unit cell 

# Andy Paul Chen, Tuesday 15 February 2022, Singapore

from Distance import * # Distance package, incl. Cartesian, Seldyn, Shmoscar
ls = longspace

def data_VacuumQD(data, vacuumSize):
    # Input: data (data file-line array) and vacuumSize (size of inter-QD spacing in angstroms)
    
    # Convert to Cartesian coordinates
    if not data_isCart(data):
        data = data_switchCart(data)
               
    # Manipulate basis vectors
    a = np.array(re.findall(r"-?\d+\.\d+", data[2].strip()))
    b = np.array(re.findall(r"-?\d+\.\d+", data[3].strip()))
    c = np.array(re.findall(r"-?\d+\.\d+", data[4].strip()))
    af = a.astype(np.float)
    bf = b.astype(np.float)
    cf = c.astype(np.float)
    sep = min(distance([0,0,0],af),distance([0,0,0],bf),distance([0,0,0],cf),distance([0,0,0],np.add(af,bf,cf)))
    print("Separation parameter: ",sep)
    
    # Cell reconstruction (big cube of side a decided by the separation parameter (min distance between lattice points)
    # and then add to that the desired gap size)
    for line in [2,3,4]:
        inp = [0.0,0.0,0.0]
        inp[line-2] = sep + vacuumSize
        data[line] = ls + "{:11f}".format(inp[0]) + ls + "{:11f}".format(inp[1]) + ls + "{:11f}".format(inp[2]) + "\n"
        
    return data
    
def VacuumQD(in_filename, out_filename, vacuumSize):
    # Take file of name in_filename, make the atoms included in the cell
    # a free-floating quantum dot separated by a vacuum of at least 
    # vaccumSize in all directions
    
    # Read input file
    data = readfile(in_filename)
    
    # Perform calculation
    dataout = data_VacuumQD(data, vacuumSize)
    
    # Write output file
    writefile(dataout, out_filename)
    