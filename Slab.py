# Slab: Create a 2D slab surrounded by vacuum from two sides, from a single unit cell 

# Andy Paul Chen, Tuesday 15 March 2022, Singapore
# Updates: The Russians have invaded Ukraine for 3 weeks now. I am going to invade church in 15 minutes.

from CQuery import * # Distance package, incl. Cartesian, Seldyn, Shmoscar
from Cartesian import *

def Slab(data, vacuumSize):
    # Input: data (data file-line array) and vacuumSize (size of inter-QD spacing in angstroms)
    # Create a 2D slab cell with a vacuum region in the c direction
    
    # Convert to Cartesian coordinates
    if not isCart(data): data = switchCart(data)
               
    # Manipulate basis vectors
    a = np.array(re.findall(r"-?\d+\.\d+", data[2].strip()))
    b = np.array(re.findall(r"-?\d+\.\d+", data[3].strip()))
    c = np.array(re.findall(r"-?\d+\.\d+", data[4].strip()))
    af = a.astype(np.float)
    bf = b.astype(np.float)
    cf = c.astype(np.float)
    
    # Norm = a cross b
    n = np.cross(af,bf)
    nu = n/np.linalg.norm(n)
    cu = cf/np.linalg.norm(cf)
    ncangle = np.arccos(np.dot(nu, cu))
    
    # Correction according to skewness of c-vector
    newnormht = np.linalg.norm(cf)*np.cos(ncangle) + vacuumSize
    ratio = newnormht/(np.linalg.norm(cf)*np.cos(ncangle))
    newc = np.multiply(cf, ratio)
    
    # Cell reconstruction (add vacuum by multiplying c-vector)
    data[4] = ls + flpr.format(newc[0]) + ls + flpr.format(newc[1]) + ls + flpr.format(newc[2]) + "\n"
        
    return data
    