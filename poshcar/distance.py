# distance: Find distance between sites in cell, in quoted Cartesian coordinates

# Andy Paul Chen, Monday, 9 August 2021, Little Italy, Cleveland, Ohio (National Day of Singapore)
# 2 Nov 2023 I can't believe I am back again

from poshcar.cartesian import * # Cartesian coordinates package

# Cite: https://doi.org/10.1039/B801115J; https://doi.org/10.1002/chem.200800987 (Bk onwards; single bond)
covalent_radius_pm = [31+5, 28, 128+7, 96+3, 84+3, 76+1, 71+1, 66+2, 57+3, 58, 166+9, 141+7, 121+4, 111+2, 107+3, 105+3, 102+4, 106+10, 203+12, 176+10, 170+7, 160+8, 153+8, 139+5, 161+8, 152+6, 150+7, 124+4, 132+4, 122+4, 122+3, 120+4, 119+4, 120+4, 120+3, 116+4, 220+9, 195+10, 190+7, 175+7, 164+6, 154+5, 147+7, 146+7, 142+7, 139+6, 145+5, 144+9, 142+5, 139+4, 139+5, 138+4, 139+3, 140+9, 244+11, 215+11, 207+8, 204+9, 203+7, 201+6, 199, 198+8, 198+6, 196+6, 194+5, 192+7, 192+7, 189+6, 190+10, 187+8, 175+10, 187+8, 170+8, 162+7, 151+7, 144+4, 141+6, 136+5, 136+6, 132+5, 145+7, 146+5, 148+4, 140+4, 150, 150, 260, 221+2, 215, 206+6, 200, 196+7, 190+1, 187+1, 180+6, 169+3, 168, 168, 165, 167, 173, 176, 161, 157, 149, 143, 141, 134, 129, 128, 121, 122, 136, 143, 162, 175, 165, 157] # add one SD

nobond = [0,2,1,1,0,0,0,0,0,2,1,1,1,0,0,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2]
# 1: Metal (M-M bond forbidden)
# 2: Noble gas (any bond forbidden)

def is_nearest(data, c, verbose = True):
    # read c (array of three floats)
    # return index of closest site to coordinates
    # display distance
    # allcoord extraction
    ns, allcoord = images(data)
    
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
                    
    if verbose: print("Minimum distance agreement: ", flpr.format(mindist), ", Index #", mindex)
    return mindex

# Andy Paul Chen, 2 November 2023
# My first coding after surviving war in Israel. They shot as us with rockets!
# Can people here ever catch a break??

def images(data):
    # Generate 27 images of an atomic site around an original unit cell (Cartesian coordinates)
    # Useful dor dealing with periodic boundary artefacts
    # Output to a unique 27-cell data structure - list of 27 lists of coordinates
    # Convert to Cartesian coordinates
    if not is_cart(data): data = switchcart(data, verbose = False)
    # Choose starting file line for coordinates according to site-disorder (selective dynamics) tag
    dcindex = 8 if is_seldyn(data) else 7
    B = basis(data) # Read lattice vectors
    
    # Read atomic coordinates
    monocoord = []
    for line in range(len(data)):
        if line > dcindex:
            dirs = np.array(re.findall(r"-?\d+\.\d+", data[line].strip()))[:3]
            dirs = dirs.astype(float) # convert from string to float
            monocoord.append(dirs)
    ns = len(monocoord) # Number of atoms in this unit cell

    # Let me write this down here real quicc. We are going to make a supercell
    # but not follow through all the processes of the supercell function. We will
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

def elemindices(data, verbose = True):
    if is_seldyn(data) and data[7].lower() == 'site-disordered structure\n': # site-disordered structures require site-specific symbols
        if verbose: print("Site-disordered structure detected!")
        dcindex = 8
        spp = []
        res = []
        occ = []
        for line in range(len(data)):
            if line > dcindex:
                spp.append(re.sub(r'[^a-zA-Z]', '', np.array(re.findall(r'\S+', data[line].strip()))[3]))
                res.append(np.array(re.findall(r'\S+', data[line].strip()))[3])
                occ.append(np.array(re.findall(r'\S+', data[line].strip()))[4])
    else:
        numatoms = np.array(re.findall('\d+', data[6].strip()))
        numatoms = numatoms.astype(int) # convert all elements to int
        elemlist = data[5].split()
        spp = list(itertools.chain.from_iterable(itertools.repeat(elemlist[i], numatoms[i]) for i in range(len(elemlist))))
        res = spp
        occ = np.ones(len(res))
    res_indexed = deepcopy(res)
    for i in range(len(res)): res_indexed[i] = res[i] + "-" + str(i+1)
    atomspp = pd.DataFrame({'Species': spp, 'Wyckoff Site': res, 'POSCAR Site': res_indexed, 'Occupancy': np.array(occ).astype(float)})
    return atomspp

def matrix_distances(data, verbose = True):
    atomspp = elemindices(data) # Header extraction
    ns, allcoord = images(data) # allcoord extraction - ns = number of atoms in cell
    
    distance_matrix = np.zeros((27,ns,ns))
    vectors_matrix = np.zeros((27,ns,ns,3))
    for i in range(ns): # i: from each atom in unit cell
        for j in range(ns): # j: to each atom in native/neighbouring cell
            for virtual in range(27): # for each virtual image
                vectors_matrix[virtual][i][j] = vector(allcoord[0][i], allcoord[virtual][j])
                distance_matrix[virtual][i][j] = distance(allcoord[0][i], allcoord[virtual][j])

    if verbose:       
        df = pd.DataFrame(distance_matrix[0], columns=list(atomspp['POSCAR Site']), index=list(atomspp['POSCAR Site']))
        print("Distances between atoms (native cell only):")
        display(df) # Display example (native cell)
    return vectors_matrix, distance_matrix  

def matrix_bonding(data, tolerance, verbose = True):
    vector_matrix, distance_matrix = matrix_distances(data, verbose = False) # Extract distance matrix, atom indices
    atomspp = elemindices(data)
    ns = len(list(atomspp['POSCAR Site']))
    covalent_radius_a = [x/100 for x in covalent_radius_pm]
    # Retrieve atomic numbers
    atomnos = np.zeros(ns).astype(int)
    for k in range(ns): atomnos[k] = periodic_table.index(list(atomspp['Species'])[k]) # extract atomic numbers
    
    # Construct bonding matrix
    # interatomic distances smaller than (radius(1) + radius(2))*(1+tolerance factor) are bonded
    bonding_matrix = np.zeros((27,ns,ns)).astype(int)
    for i in range(ns): # i: from atom
        for j in range(ns): # j: to atom
            ei = atomnos[i] # atomic number of species i
            ej = atomnos[j] # atomic number of species j
            for virtual in range(27): # for each virtual image
                if distance_matrix[virtual][i][j] <= (covalent_radius_a[ei]+covalent_radius_a[ej])*(tolerance+1):
                    if (i != j) and not (nobond[ei]==1 and nobond[ej]==1) and not(nobond[ei]==2 or nobond[ej]==2):
                        bonding_matrix[virtual][i][j] = 1
    
    if verbose:
        df = pd.DataFrame(bonding_matrix[0], columns=list(atomspp['POSCAR Site']), index=list(atomspp['POSCAR Site']))
        print("Bonding between atoms (native cell only):")
        display(df) # Display example (native cell)
    return bonding_matrix

def matrix_bonding_average(data, mode, tolerance, bme_correlated = 'amaiwana', verbose = True):
    # Average bonding matrix, can also be used for site-disordered stuff
    # mode: string input, first letter for [s]ite classification or [e]lement species classification
    # Edit of 10 Sep 2024: deviations from the mean are recorded in a separate matrix to output
    # bme_correlated specified if another averaged matrix needs to be used (this is in the case of correlated disorder)

    atomspp = elemindices(data, verbose) # What atoms are in here
    weights = atomspp['Occupancy']
    bms = matrix_bonding(data, tolerance, verbose = False).sum(0)

    # Mode 
    classification = mode[0].upper()
    if classification == 'E': arrayclass = np.array(atomspp['Species'])
    elif classification == 'S': arrayclass = np.array(atomspp['Wyckoff Site'])
    else: 
        print("Invalid mode -- please enter 'site' or 'element'! Defaulting to element species classification...\n")
        arrayclass = np.array(atomspp['Species'])
    
    runiq, runiq_index, runiq_count = np.unique(arrayclass, return_counts=True, return_index=True)
    runiq = runiq[np.argsort(runiq_index)]
    runiq_count = runiq_count[np.argsort(runiq_index)] # preserve order
    us = len(runiq) # How many unique species in cell
    # target bme
    bme = np.zeros([us, us])
    if type(bme_correlated) != str: bmet = bme_correlated
    
    res_partition = np.resize(np.cumsum(runiq_count),runiq_count.size -1) # partition plan of array
    bm_split_by_rows = np.split(bms,res_partition,axis=0)
    weights_split = np.split(weights,res_partition)
    stackerv_bma = []
    stackerv_bmc = []
    for iu in range(us): # rows split
        bm_split_by_columns = np.split(bm_split_by_rows[iu],res_partition,axis=1)
        stackerh_bma = []
        stackerh_bmc = []
        for ju in range(us): # iterate by column
            summed_column = np.dot(bm_split_by_columns[ju],weights_split[ju])
            bme[iu][ju] = np.average(summed_column, weights = weights_split[iu])
            summed_column = summed_column.reshape(-1,1)
            piece = summed_column.copy() # this piece goes into bma
            if type(bme_correlated) == str: piece.fill(bme[iu][ju]) # make a matrix full of the same number
            else: piece.fill(bmet[iu][ju]) # use target to make bma
            piece = piece.reshape(-1,1)
            stackerh_bma.append(piece) # stack horizontal bma
            stackerh_bmc.append(summed_column)
        slab_bma = np.hstack(stackerh_bma)
        slab_bmc = np.hstack(stackerh_bmc)
        stackerv_bma.append(slab_bma) # stack vertical
        stackerv_bmc.append(slab_bmc) # stack vertical
    # bme is made by this point

    bma = np.vstack(stackerv_bma) # bma is made here
    bmc = np.vstack(stackerv_bmc) # bmc is made here
    
    # Display and return bme
    if verbose:
        # Display and return bme
        df = pd.DataFrame(bme, columns=runiq, index=runiq)
        print("bme - average local coordination by environment:")
        display(df) # Display example (native cell)

        # Display and return bma
        df1 = pd.DataFrame(bma, columns=runiq, index=atomspp['POSCAR Site'])
        print("bma - average local coordination by environment, expanded")
        display(df1) # Display example (native cell)


        # Display and return bmc
        df2 = pd.DataFrame(bmc, columns=runiq, index=atomspp['POSCAR Site'])
        print("bme - average local coordination by POSCAR site:")
        display(df2) # Display example (native cell)
    
    # Calculate "unaverageness" aka Pauling-5 penalty 
    unaverageness = np.linalg.norm(bma-bmc, 'fro')
    print("Unaverageness: ", unaverageness)

    return runiq, bme, unaverageness

def crashtest(data, tolerance, verbose = True):
    # Determine if an atom pair is too close to one another
    # tolerance is a percentage value (0.0-1.0)
    # 26 Jul 2024: returns crash matrix instead of bool output

    tolerance = min(abs(tolerance),1.0)
    if verbose: print ("Tolerance: ", tolerance)
    accept = True
    
    vector_matrix, distance_matrix = matrix_distances(data, verbose) # Extract distance matrix, atom indices
    atomspp = elemindices(data, verbose)
    ns = len(list(atomspp['POSCAR Site']))
    covalent_radius_a = [x/100 for x in covalent_radius_pm]
    # Retrieve atomic numbers
    atomnos = np.zeros(ns).astype(int)
    for k in range(ns): atomnos[k] = periodic_table.index(list(atomspp['Species'])[k]) # extract atomic numbers
    
    # Construct bonding matrix
    # interatomic distances smaller than (radius(1) + radius(2))*(1+tolerance factor) are bonded
    crashing_matrix = np.zeros((ns,ns)).astype(int)
    
    for i in range(ns): # i: from atom
        for j in range(ns): # j: to atom
            ei = atomnos[i] # atomic number of species i
            ej = atomnos[j] # atomic number of species j
            for virtual in range(27): # for each virtual image
                if distance_matrix[virtual][i][j] <= (covalent_radius_a[ei]+covalent_radius_a[ej])*(1-tolerance):
                    if (i != j):
                        crashing_matrix[i][j] = 1
                        accept = False
                    
    df = pd.DataFrame(crashing_matrix, columns=list(atomspp['POSCAR Site']), index=list(atomspp['POSCAR Site']))
    if verbose:
        if accept == False: 
            print("Crashing between atoms:")
            display(df) # Display example
        else: print("No crashing between atoms!")
    return crashing_matrix