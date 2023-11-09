# VacuumQD: Create a QD surrounded by vacuum from a single unit cell 

# Andy Paul Chen, Tuesday 15 February 2022, Singapore
# (Update of 10 May 2022): Termination control (same termination, all round!) in VacuumQDPerov

from Distance import * # Distance package, incl. Cartesian, Seldyn, Shmoscar
from Perovskite import *
from Clip import *
from Supercell import *

def VacuumQD(data, vacuumSize):
    # Input: data (data file-line array) and vacuumSize (size of inter-QD spacing in angstroms)
    
    # Convert to Cartesian coordinates
    if not isCart(data): data = switchCart(data)
               
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
    
def VacuumQDPerov(a, perovskite, supercell, vacuumSize, termination):
    # returns a perovskite quantum dot
    # of formula ABX_3 (perovskite = ["A", "B", "X"])
    # and cells reduplicated (supercell=[n1,n2,n3]) times in the [a,b,c] directions
    # with (vacuumSize) spacing between cells
    # Termination ("A" or "B") specifies whether AO or BO2 termination used
    # If neither is specified, the termination is not adjusted
    
    # Generate base perovskite unit cell
    data = perovABX3(a, perovskite[0], perovskite[1], perovskite[2])
    
    # Generate Supercell
    data = Supercell(data, supercell[0], supercell[1], supercell[2])
    
    # Insert vacuum spacing
    data = VacuumQD(data, vacuumSize)
    
    # Enforce termination
    marg = a/4 # defines margin for clipping
    if termination.upper() == "A":
        # clip at far edge
        cutoff = [(a*supercell[0])-(3*marg), (a*supercell[1])-(3*marg), (a*supercell[2])-(3*marg)]
        data = Clip(data, 1, cutoff)
    elif termination.upper() == "B":
        # clip at near edge
        cutoff = [marg, marg, marg]
        data = Clip(data, -1, cutoff)
    else:
        print("WARNING: No valid input specified for termination (QD has mixed termination)! \n")
    
    return data
    