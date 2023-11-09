# Matrix: Provides a matrix of all interatomic distances given a unit cell

# Andy Paul Chen, 2 November 2023
# My first coding after surviving war in Israel. They shot as us with rockets!
# Can people here ever catch a break??

from Distance import *

def Matrix_Distances(data):
    # Convert to Cartesian coordinates, check for selective dynamics
    if not isCart(data): data = switchCart(data)
    # Cancel selective dynamics if active
    if isSeldyn(data): SeldynSwitch(data)
    dcindex = 7
    
    # Read lattice vectors
    m1 = re.findall(r"-?\d+\.\d+", data[1].strip())
    mu = float(m1[0])
    a = re.findall(r"-?\d+\.\d+", data[2].strip())
    a = [float(x)*mu for x in a]
    b = re.findall(r"-?\d+\.\d+", data[3].strip())
    b = [float(x)*mu for x in b]
    c = re.findall(r"-?\d+\.\d+", data[4].strip())
    c = [float(x)*mu for x in c]
    
    # Header extraction
    numatoms = np.array(re.findall('\d+', data[6].strip()))
    numatoms = numatoms.astype(int) # convert all elements to int
    elemlist = data[5].split()
    res = list(itertools.chain.from_iterable(itertools.repeat(elemlist[i], numatoms[i]) for i in range(len(elemlist))))
    res_indexed = deepcopy(res)
    for i in range(len(res)): res_indexed[i] = res[i] + "-" + str(i+1)
    print("List of atoms: ", res_indexed)
    
    # Read atomic coordinates
    monocoord = []
    for line in range(len(data)):
        if line > dcindex:
            dirs = re.findall(r"-?\d+\.\d+", data[line].strip())
            dirs = [float(x) for x in dirs] # convert from string to float
            monocoord.append(dirs)
    nosites = len(monocoord) # Number of atoms in this unit cell
    
    # Let me write this down here real quicc. We are going to make a supercell
    # but not follow through all the processes of the Supercell function. We will
    # extend the list of atoms to the 27 surrounding unit cells only. Closest
    # neighbours which are images on the list will be mapped to the original ion.
    
    # Make 26 shadow-coordinate cells
    allcoord = [deepcopy(monocoord)]
    for it in range(26): allcoord += [deepcopy(monocoord)]
    
    # Move 26 shadow-coordinate cells to their respective positions
    for i in [9,10,11,12,13,14,15,16,17]:
        for j in range(nosites):
            allcoord[i][j] = list(map(add, allcoord[i][j], a))
    for i in [18,19,20,21,22,23,24,25,26]:
        for j in range(nosites):
            allcoord[i][j] = list(map(add, allcoord[i][j], [-x for x in a]))
    for i in [3,4,5,12,13,14,21,22,23]:
        for j in range(nosites):
            allcoord[i][j] = list(map(add, allcoord[i][j], b))
    for i in [6,7,8,15,16,17,24,25,26]:
        for j in range(nosites):
            allcoord[i][j] = list(map(add, allcoord[i][j], [-x for x in b]))
    for i in [1,4,7,10,13,16,19,22,25]:
        for j in range(nosites):
            allcoord[i][j] = list(map(add, allcoord[i][j], c))
    for i in [2,5,8,11,14,17,20,23,26]:
        for j in range(nosites):
            allcoord[i][j] = list(map(add, allcoord[i][j], [-x for x in c]))
    
    # Create distances matrix
    distance_matrix = np.zeros((nosites,nosites))
    for i in range(nosites): # i: from atom
        for j in range(nosites): # j: to atom
            for virtual in range(27):
                distance_to_virtual = distance(allcoord[0][i], allcoord[virtual][j])
                if virtual==0: min_dist = distance_to_virtual
                if (distance_to_virtual <= min_dist) or (min_dist ==0): min_dist = distance_to_virtual
                distance_matrix[i][j] = min_dist 
    df = pd.DataFrame(distance_matrix, columns=res_indexed, index=res_indexed)
    return df