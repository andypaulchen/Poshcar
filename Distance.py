# Distance: Find distance between sites in cell, in quoted Cartesian coordinates

# Andy Paul Chen, Monday, 9 August 2021, Little Italy, Cleveland, Ohio (National Day of Singapore)
# 2 Nov 2023 I can't believe I am back again

from Cartesian import * # Cartesian coordinates package
from CQuery import * # coordinates query package

covalent_radius_pm = [25, 0, 145, 105, 85, 70, 65, 60, 50, 0, 180, 150, 125, 110, 100, 100, 100, 0, 220, 180, 160, 140, 135, 140, 140, 140, 135, 135, 135, 135, 130, 125, 115, 115, 115, 0, 235, 200, 180, 155, 145, 145, 135, 130, 135, 140, 160, 155, 155, 145, 145, 140, 140, 0, 260, 215, 195, 185, 185, 185, 185, 185, 185, 180, 175, 175, 175, 175, 175, 175, 175, 155, 145, 135, 135, 130, 135, 135, 135, 150, 190, 180, 160, 190, 0, 0, 0, 215, 195, 180, 180, 175, 175, 175, 175, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] 

def distance(c1, c2):
    # Read Cartesian coordinates 1 and 2, derive distance (in Angstroms)
    # c1, c2 are np. arrays of three float values
    a = np.array(c1)
    b = np.array(c2)
    return np.linalg.norm(a-b)

def vector(c1, c2):
    a = np.array(c1)
    b = np.array(c2)
    return b-a

def isNearest(data, c):
    # read c (array of three floats)
    # return index of closest site to coordinates
    # display distance
    # allcoord extraction
    ns, allcoord = Images(data)
    
    # analysis
    for i in range(ns):
        for virtual in range(27):
            dist = distance(np.array(c), allcoord[virtual][i])
            if i==0 and virtual==0:
                mindex = 1
                mindist = dist
            else:
                if dist < mindist:
                    mindex = i+1
                    mindist = dist
                    
    print("Minimum distance agreement: ", flpr.format(mindist), ", Index #", mindex)
    return mindex

# Andy Paul Chen, 2 November 2023
# My first coding after surviving war in Israel. They shot as us with rockets!
# Can people here ever catch a break??

def Images(data):
    # Generate 27 images of an atomic site around an original unit cell (Cartesian coordinates)
    # Useful dor dealing with periodic boundary artefacts
    # Output to a unique 27-cell data structure - list of 27 lists of coordinates
    # Convert to Cartesian coordinates
    if not isCart(data): data = switchCart(data)
    # Cancel selective dynamics if active
    if isSeldyn(data): data = SeldynSwitch(data)
    dcindex = 7
    B = Basis(data) # Read lattice vectors
    
    # Read atomic coordinates
    monocoord = []
    for line in range(len(data)):
        if line > dcindex:
            dirs = np.array(re.findall(r"-?\d+\.\d+", data[line].strip()))
            dirs = dirs.astype(float) # convert from string to float
            monocoord.append(dirs)
    ns = len(monocoord) # Number of atoms in this unit cell

    # Let me write this down here real quicc. We are going to make a supercell
    # but not follow through all the processes of the Supercell function. We will
    # extend the list of atoms to the 27 surrounding unit cells only. Closest
    # neighbours which are images on the list will be mapped to the original ion.
    
    # Make 26 shadow-coordinate cells
    allcoord = [deepcopy(monocoord)]
    for it in range(26): allcoord += [deepcopy(monocoord)]
    
    # Move 26 shadow-coordinate cells to their respective positions
    for i in [9,10,11,12,13,14,15,16,17]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] + B[0] # +a
    for i in [18,19,20,21,22,23,24,25,26]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] - B[0] # -a
    for i in [3,4,5,12,13,14,21,22,23]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] + B[1] # +b
    for i in [6,7,8,15,16,17,24,25,26]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] - B[1] # -b
    for i in [1,4,7,10,13,16,19,22,25]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] + B[2] # +c
    for i in [2,5,8,11,14,17,20,23,26]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] - B[2] # -c
            
    return ns, allcoord

def ElemIndices(data):
    numatoms = np.array(re.findall('\d+', data[6].strip()))
    numatoms = numatoms.astype(int) # convert all elements to int
    elemlist = data[5].split()
    res = list(itertools.chain.from_iterable(itertools.repeat(elemlist[i], numatoms[i]) for i in range(len(elemlist))))
    res_indexed = deepcopy(res)
    for i in range(len(res)): res_indexed[i] = res[i] + "-" + str(i+1)
    print("List of atoms: ", res_indexed)
    return res, res_indexed

def Matrix_Distances(data):
    # Header extraction
    res, res_indexed = ElemIndices(data)
    # allcoord extraction
    ns, allcoord = Images(data)
    
    # Create distances matrix
    distance_matrix = np.zeros((ns,ns))
    vectors_matrix = []
    for i in range(ns): # i: from atom
        vectors_list = []
        for j in range(ns): # j: to atom
            #print(res_indexed[i], res_indexed[j], i, j)
            for virtual in range(27):
                vector_to_virtual = vector(allcoord[0][i], allcoord[virtual][j])
                distance_to_virtual = distance(allcoord[0][i], allcoord[virtual][j])
                #print("\t Image ID: ", virtual, vector_to_virtual, distance_to_virtual)
                if virtual==0: 
                    min_dist = distance_to_virtual
                    min_vec = vector_to_virtual[:]
                    #print("\t\t min_dist: ", min_dist, "/ min_vec : ", min_vec)
                if (distance_to_virtual <= min_dist) or (min_dist == 0): 
                    min_dist = distance_to_virtual
                    min_vec = vector_to_virtual[:]
                    #print("\t\t min_dist: ", min_dist, "/ min_vec : ", min_vec)
            distance_matrix[i][j] = min_dist 
            vectors_list.append(min_vec)
        vectors_matrix.append(vectors_list)
                
    df = pd.DataFrame(distance_matrix, columns=res_indexed, index=res_indexed)
    return vectors_matrix, distance_matrix, df

def Matrix_Bonding(data, tolerance):
    # Extract distance matrix, atom indices
    vector_matrix, distance_matrix, df = Matrix_Distances(data)
    res, res_indexed = ElemIndices(data)
    ns = len(res_indexed)
    covalent_radius_a = [x/100 for x in covalent_radius_pm]
    
    # Retrieve atomic numbers
    atomnos = np.zeros(ns).astype(int)
    for k in range(ns):
        atomnos[k] = Periodic_Table.index(res[k])
    
    # Construct bonding matrix
    # interatomic distances smaller than (radius(1) + radius(2))*(1+tolerance factor) are bonded
    bonding_matrix = np.zeros((ns,ns)).astype(int)
    for i in range(ns): # i: from atom
        for j in range(ns): # j: to atom
            if distance_matrix[i][j] <= (covalent_radius_a[atomnos[i]]+covalent_radius_a[atomnos[j]])*(tolerance+1):
                bonding_matrix[i][j] = 1
    df = pd.DataFrame(bonding_matrix, columns=res_indexed, index=res_indexed)
    return bonding_matrix, df