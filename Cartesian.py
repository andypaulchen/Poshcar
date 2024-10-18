# Cartesian: Convert Direct to Cartesian coordinates, and vice versa (?)

# Andy Paul Chen, Tuesday 21 April 2020, Cleveland Heights, Ohio

# 18 February 2022: updated to account for line 1 multiplier of lattice vectors

from Seldyn import * # selective dynamics package

def isCart(data):
    dcindex = 8 if isSeldyn(data) else 7 # Index of Direct/Cartesian line
    if data[dcindex][0].upper() == 'C' or data[dcindex][0].upper() == 'K':
        return True
        
def switchCart(data, verbose = True):
    # Switches coordinates between Cartesian and Direct!
    dcindex = 8 if isSeldyn(data) else 7 # Index of Direct/Cartesian line
    num = re.findall(r'\d+', data[atom_number_index].strip())
    snum = sum(list(map(int,num)))
    data = data[:(snum + dcindex + 1)] # whittle the tail
    B = Basis(data) # Read lattice vectors
    
    # Announce the conversion operation
    if verbose:
        if isCart(data): print("Converting: Cartesian > Fractional")
        else: print("Converting: Fractional > Cartesian")
        
    # Conversion operation
    for line in range(len(data)):
        if line > dcindex:
            F = np.array(re.findall(r"-?\d+\.\d+", data[line].strip()))
            # Read seldyn flags for later
            if isSeldyn(data): flags = re.findall(r"\w", data[line].strip())[-3:]
            # Convert the line - MATRIX MULTIPLICATION
            if isCart(data): outv = np.matmul(np.linalg.inv(B.transpose()), F.astype(float)) % 1
            else: outv = np.matmul(B.transpose(), F.astype(float))
            # write line
            data[line] = ls + flpr.format(outv[0]) + ls + flpr.format(outv[1]) + ls + flpr.format(outv[2])
            # Add seldyn flags back in
            if isSeldyn(data): data[line] += ls + flags[0] + " " + flags[1] + " " + flags[2]
            data[line] += "\n"
    
    # Switch the Direct/Cartesian line
    data[dcindex] = "Direct\n" if isCart(data) else "Cartesian\n"
    return data

def Translate(data, vector):
    # Convert to Cartesian coordinates
    if not isCart(data): data = switchCart(data)
    if isSeldyn(data): data = SeldynSwitch(data)
    vector = np.array(vector)
    vector = vector.astype(float)
    dcindex = 8 if isSeldyn(data) else 7 # Index of Direct/Cartesian line
    for line in range(len(data)):
        if line > dcindex:
            coords = np.array(re.findall(r"-?\d+\.\d+", data[line].strip()))
            coords = coords.astype(float)
            newc = coords + vector
            data[line] = ls + str(newc[0]) + ls + str(newc[1]) + ls + str(newc[2]) + "\n"
    return data