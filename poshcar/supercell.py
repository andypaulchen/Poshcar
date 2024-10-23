# supercell: Generates a repeated supercell from a smaller cell

# Andy Paul Chen, Monday 21 February 2022, Singapore
# 9 May 2024, added the line on "tailmatter" so supercells maintain seldyn / disorder info

from poshcar.cartesian import * 

def supercell(data, AA, BB, CC):
    if (type(AA) is int) and (type(BB) is int) and (type(CC) is int):
        # Convert to Cartesian coordinates
        if not is_cart(data): data = switchcart(data)
        dcindex = 8 if is_seldyn(data) else 7 # Index of Direct/Cartesian line
            
        # Manipulate basis vectors
        B = basis(data) # Read lattice vectors
        # Multiply by multipliers
        aff = np.multiply(B[0],AA)
        bff = np.multiply(B[1],BB)
        cff = np.multiply(B[2],CC)
        # Rewrite lines
        data[2] = ls + flpr.format(aff[0]) + ls + flpr.format(aff[1]) + ls + flpr.format(aff[2]) + "\n"
        data[3] = ls + flpr.format(bff[0]) + ls + flpr.format(bff[1]) + ls + flpr.format(bff[2]) + "\n"
        data[4] = ls + flpr.format(cff[0]) + ls + flpr.format(cff[1]) + ls + flpr.format(cff[2]) + "\n"
        
        # Manipulate atom species numbers
        numatoms = np.array(re.findall('\d+', data[6].strip()))
        numatoms = [int(i) for i in numatoms] # convert all elements to int
        numnew = np.multiply(numatoms,AA*BB*CC)
        linenew = ""
        for num in numnew:
            linenew = linenew + ls + str(num)
        linenew = linenew + "\n"
        data[6] = linenew
        
        # Split POSCAR into header and element-specific segments
        header = data[:(dcindex+1)] # entire file except atomic coordinate lines
        rem = data[(dcindex+1):] # the atomic coordinate lines
        splitlists = []
        for count in numatoms:
            splitlists.append(rem[:count])
            rem = rem[count:]
        
        # Extend
        for elem in list(range(len(numatoms))):
            sublist = splitlists[elem] # list containing all coordinates of particular elem
            sub2list = sublist.copy()
            tailmatter = sublist.copy()
            # convert to float
            for oo in list(range(len(sublist))):
                hll = re.findall(r'\S+', sublist[oo].strip())
                tailmatter[oo] = ' '.join(hll[3:]) # material beyond the coordinates
                sub2list[oo] = [float(hll[0]),float(hll[1]),float(hll[2])]
            newlist = []
            
            # replicate atoms by AA*BB*CC times
            for ii in list(range(AA)):
                for jj in list(range(BB)):
                    for kk in list(range(CC)):
                        translate = np.add(np.add(np.multiply(ii,B[0]),np.multiply(jj,B[1])),np.multiply(kk,B[2]))
                        replicate = sub2list.copy()
                        for pp in list(range(len(replicate))):
                            nc = np.add(sub2list[pp], translate) # new coordinates
                            replicate[pp] = ls + flpr.format(nc[0]) + ls + flpr.format(nc[1]) + ls + flpr.format(nc[2]) + ls + tailmatter[pp] + "\n"
                        newlist.extend(replicate)

            # Integrate head and body matter
            header += newlist
        return header
    else:
        print("Please input integers as arguments!!\n")