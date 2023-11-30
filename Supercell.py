# Supercell: Generates a repeated supercell from a smaller cell

# Andy Paul Chen, Monday 21 February 2022, Singapore

from Cartesian import * # SHMOSCAR core package

def Supercell(data, AA, BB, CC):
    if (type(AA) is int) and (type(BB) is int) and (type(CC) is int):
        # Convert to Cartesian coordinates
        if not isCart(data): data = switchCart(data)
        dcindex = 8 if isSeldyn(data) else 7 # Index of Direct/Cartesian line
            
        # Manipulate basis vectors
        B = Basis(data) # Read lattice vectors
        # Multiply by multipliers
        aff = np.multiply(B[0],AA)
        bff = np.multiply(B[1],BB)
        cff = np.multiply(B[2],CC)
        # Rewrite lines
        data[2] = ls + "{:11f}".format(aff[0]) + ls + "{:11f}".format(aff[1]) + ls + "{:11f}".format(aff[2]) + "\n"
        data[3] = ls + "{:11f}".format(bff[0]) + ls + "{:11f}".format(bff[1]) + ls + "{:11f}".format(bff[2]) + "\n"
        data[4] = ls + "{:11f}".format(cff[0]) + ls + "{:11f}".format(cff[1]) + ls + "{:11f}".format(cff[2]) + "\n"
        
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
        header = data[:(dcindex+1)]
        rem = data[(dcindex+1):]
        splitlists = []
        for count in numatoms:
            splitlists.append(rem[:count])
            rem = rem[count:]
        
        # Extend
        for elem in list(range(len(numatoms))):
            sublist = splitlists[elem]
            sub2list = sublist.copy()
            # convert to float
            for oo in list(range(len(sublist))):
                hll = re.findall(r"-?\d+\.\d+", sublist[oo].strip())
                sub2list[oo] = [float(hll[0]),float(hll[1]),float(hll[2])]
            newlist = []
            
            # replicate atoms by AA*BB*CC times
            for ii in list(range(AA)):
                for jj in list(range(BB)):
                    for kk in list(range(CC)):
                        translate = np.add(np.add(np.multiply(ii,B[0]),np.multiply(jj,B[1])),np.multiply(kk,B[2]))
                        replicate = sub2list.copy()
                        for pp in list(range(len(replicate))):
                            replicate[pp] = np.add(sub2list[pp], translate)
                        newlist.extend(replicate)
                        
            for qq in list(range(len(newlist))):
                ph = newlist[qq]
                newlist[qq] = ls + "{:11f}".format(ph[0]) + ls + "{:11f}".format(ph[1]) + ls + "{:11f}".format(ph[2]) + "\n"
            header = header + newlist
        return header
    else:
        print("Please input integers as arguments!!\n")