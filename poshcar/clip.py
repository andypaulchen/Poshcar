# clip: Deletes atoms beyond/before a certain coordinate value
# (uses atomsub to delete atoms)

# Andy Paul Chen, Monday 9 May 2022

from poshcar.atomsub import *
from poshcar.cartesian import *

def clip(data, greaterorsmallerthan, cutoff):    
    # In the case of Selective Dynamics, the first coordinates are on line 8
    # Otherwise, it's 7
    dcindex = 8 if is_seldyn(data) else 7
    
    # Change everything into Cartesian
    if not is_cart(data): data = switchcart(data, verbose = False)
    
    # Find all atoms of index [#] to be deleted into deletelist
    trackindex = 0
    deletelist = [0]
    for line in range(len(data)):
        if line > dcindex:
            trackindex = trackindex + 1
            dirs = re.findall(r"-?\d+\.\d+", data[line].strip())
            if greaterorsmallerthan > 0:
                if (float(dirs[0]) > cutoff[0]) or (float(dirs[1]) > cutoff[1]) or (float(dirs[2]) > cutoff[2]):
                    deletelist.append(trackindex)
            else:
                if (float(dirs[0]) < cutoff[0]) or (float(dirs[1]) < cutoff[1]) or (float(dirs[2]) < cutoff[2]):
                    deletelist.append(trackindex)
    # Remove initializing element of deletelist
    deletelist.remove(0)
    # You want deletelist in descending order, давай?
    deletelist.sort(reverse = True)
    
    # Remove the atoms
    data = atomsub(data, "vac", deletelist, verbose = False)
    return data