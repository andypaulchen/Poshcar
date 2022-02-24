# Supercell: Generates a repeated supercell from a smaller cell

# Andy Paul Chen, Monday 21 February 2022, Singapore

from Cartesian import * # SHMOSCAR core package
ls = longspace
lls = ls + "   "

def data_Supercell(data, AA, BB, CC):
    # Convert to Cartesian coordinates
    if not data_isCart(data):
        data = data_switchCart(data)
    dcindex = 8 if data_isSeldyn(data) else 7 # Index of Direct/Cartesian line
        
    # Manipulate basis vectors
    a = np.array(re.findall(r"-?\d+\.\d+", data[2].strip()))
    b = np.array(re.findall(r"-?\d+\.\d+", data[3].strip()))
    c = np.array(re.findall(r"-?\d+\.\d+", data[4].strip()))
    # Convert data to float
    af = a.astype(np.float)
    bf = b.astype(np.float)
    cf = c.astype(np.float)
    # Multiply by multipliers
    aff = np.multiply(af,AA)
    bff = np.multiply(bf,BB)
    cff = np.multiply(cf,CC)
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
                    translate = np.add(np.add(np.multiply(ii,af),np.multiply(jj,bf)),np.multiply(kk,cf))
                    replicate = sub2list.copy()
                    for pp in list(range(len(replicate))):
                        replicate[pp] = np.add(sub2list[pp], translate)
                    newlist.extend(replicate)
                    
        for qq in list(range(len(newlist))):
            ph = newlist[qq]
            newlist[qq] = ls + "{:11f}".format(ph[0]) + ls + "{:11f}".format(ph[1]) + ls + "{:11f}".format(ph[2]) + "\n"
        header = header + newlist
    
    return header
    
    
def Supercell(in_filename, out_filename, AA, BB, CC):
    # Takes POSCAR from in_filename, build a supercell which is the original cell
    # duplicated by AA times in the a direction, and so on.
    
    # Read input file
    data = readfile(in_filename)
    
    if (type(AA) is int) and (type(BB) is int) and (type(CC) is int):
        data = data_Supercell(data, AA, BB, CC)
    else:
        print("Please input integers as arguments!!\n")
    
    # Write output file
    writefile(data, out_filename)