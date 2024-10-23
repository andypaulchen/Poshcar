# slab: Create a 2D slab surrounded by vacuum from two sides, from a single unit cell 

# Andy Paul Chen, Tuesday 15 March 2022, Singapore
# Updates: The Russians have invaded Ukraine for 3 weeks now. I am going to invade church in 15 minutes.

from poshcar.cquery import * # Distance package, incl. Cartesian, Seldyn, Shmoscar
from poshcar.cartesian import *

def slab(data, vacuumSize):
    # Input: data (data file-line array) and vacuumSize (size of inter-QD spacing in angstroms)
    # Create a 2D slab cell with a vacuum region in the c direction
    
    # Convert to Cartesian coordinates
    if not is_cart(data): data = switchcart(data, verbose = False)
               
    # Manipulate basis vectors
    a = np.array(re.findall(r"-?\d+\.\d+", data[2].strip()))
    b = np.array(re.findall(r"-?\d+\.\d+", data[3].strip()))
    c = np.array(re.findall(r"-?\d+\.\d+", data[4].strip()))
    af = a.astype(float)
    bf = b.astype(float)
    cf = c.astype(float)
    
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

def yoink(data, n, yoinkdist):
    # moves a group of atoms (index >= n) in the +c direction in slab model

    # Convert to Cartesian coordinates
    if not is_cart(data): data = switchcart(data, verbose = False)
    dcindex = 8 if is_seldyn(data) else 7 # Index of Direct/Cartesian line

    for line in range(len(data)):
        if line >= dcindex+n:
            F = np.array(re.findall(r"-?\d+\.\d+", data[line].strip()))
            F = F.astype(float)
            # Read seldyn flags for later
            if is_seldyn(data): flags = re.findall(r"\w", data[line].strip())[-3:]
            # Write line
            F[2] = F[2] + yoinkdist
            newline = ls + flpr.format(F[0]) + ls + flpr.format(F[1]) + ls + flpr.format(F[2])
            # Add seldyn flags back in
            if is_seldyn(data): newline += ls + flags + "\n"
            else: newline += "\n"
            data[line] = newline

    return data