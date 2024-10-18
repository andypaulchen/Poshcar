# Additive.py: Add an atom to the cell

# 1 December 2023: WHY DIDN'T I DO THIS EARLIER? I had a function to remove an atom (create vacancy)
# but not to add an atom. Silly me. Maybe I can make interstitial defects here also.

from Cartesian import *
from CQuery import *
from AtomSub import *

def Tidy(data):
    # Makes file look tidier by setting all floats to same width
    for i in range(2,5):
        bv = np.array(re.findall(r"-?\d+\.\d+", data[i].strip()))
        bv = bv.astype(float)
        data[i] = ls + flpr.format(bv[0]) + ls + flpr.format(bv[1]) + ls + flpr.format(bv[2]) + "\n"
    
    dcindex = 8 if isSeldyn(data) else 7 # Index of Direct/Cartesian line
    for line in range(len(data)):
        if line > dcindex:
            F = np.array(re.findall(r"-?\d+\.\d+", data[line].strip()))
            F = F.astype(float)
            # Read seldyn flags for later
            if isSeldyn(data): flags = re.sub(r'^\s*[\d.\s]+', '', data[line])
            # Write line
            newline = ls + flpr.format(F[0]) + ls + flpr.format(F[1]) + ls + flpr.format(F[2])
            # Add seldyn flags back in
            if isSeldyn(data): newline += ls + flags
            else: newline += "\n"
            data[line] = newline
    
    return data

def AddAtom(data, elem, coords, verbose = True):
    addindex = 8 if isSeldyn(data) else 7 # Index of Direct/Cartesian line
    # Direct or Cartesian? Coordinates just follow the preset regardless
    if verbose:
        if isCart(data): 
            print("Adding "+elem+" atom, Cartesian coordinates ["+str(coords[0])+", "+str(coords[1])+", "+ str(coords[2])+"]\n")
        else: print("Adding "+elem+" atom, Fractional coordinates ["+str(coords[0])+", "+str(coords[1])+", "+ str(coords[2])+"]\n")
    
    # Read string of element symbols and number of atoms per element
    elem_list = re.findall(r'\w+', data[atom_name_index].strip())
    temp = re.findall(r'\d+', data[atom_number_index].strip())
    num_list = list(map(int, temp))
    
    # Append atom in element lists
    data[atom_number_index] = ""
    data[atom_name_index] = ""
    if elem in elem_list: 
        num_list[elem_list.index(elem)] += 1
        newatomindex = sum(num_list[:(elem_list.index(elem)+1)])
    else: 
        elem_list += [elem]
        num_list += [1]
        newatomindex = sum(num_list)
    
    for m in range(len(num_list)):
        data[atom_name_index] += ls + elem_list[m]
        data[atom_number_index] += ls + str(num_list[m])
    data[atom_name_index] += "\n" 
    data[atom_number_index] += "\n"  
    
    newline = ls + flpr.format(coords[0]) + ls + flpr.format(coords[1]) + ls + flpr.format(coords[2])
    newline += "\n" # Don't add seldyn flags anymore
    data.insert(addindex + newatomindex, newline)
    
    data = Tidy(data)
    if verbose: printvaspdata(data)
    return data

def RemoveAtom(data, indices, verbose = True):
    data = AtomSub(data,"vac",indices, verbose)
    # Read string of element symbols and number of atoms per element
    elem_list = re.findall(r'\w+', data[atom_name_index].strip())
    temp = re.findall(r'\d+', data[atom_number_index].strip())
    num_list = list(map(int, temp))
    data[atom_name_index] = data[atom_number_index] = ""
    for i in range(len(num_list)):
        if num_list[i] != 0:
            data[atom_name_index] += ls + elem_list[i]
            data[atom_number_index] += ls + str(num_list[i])
    data[atom_name_index] += "\n" 
    data[atom_number_index] += "\n" 
    
    data = Tidy(data)
    if verbose: printvaspdata(data)
    return data

def Graft(data_base, point_base, data_graft, point_graft):
    # Grafts a set of atoms ("graft" e.g. a molecule) onto a basic set ("base", e.g. slab or QD)
    # This is intended to be implemented to add ligands to a surface
    # Change both cells to CARTESIAN coordinates
    if not isCart(data_base): data_base = switchCart(data_base, verbose = False)
    if not isCart(data_graft): data_graft = switchCart(data_graft, verbose = False)
    
    # Handle translation
    data_graft = Translate(data_graft, np.array(point_base)-np.array(point_graft))
    
    dcindex_graft = 8 if isSeldyn(data_graft) else 7 # Index of Direct/Cartesian line
    for line in range(len(data_graft)):
        if line > dcindex_graft:
            elem, coord = CQuery(data_graft, line-dcindex_graft, verbose = False)
            to_add = np.array(re.findall(r"-?\d+\.\d+", data_graft[line].strip()))
            to_add = to_add.astype(float)
            data_base = AddAtom(data_base, elem, to_add, verbose = False)
    
    #Preview of grafted cell
    #writefile(data_base, "_operations/___graft.vasp")
    #view(read("_operations/___graft.vasp"), viewer='ngl')
    #user = input("Proceed with graft? (type [y]es to approve): ")
    #if user.lower()[0] == 'y': return data_base
    #else: print("Try with other parameters!\n")
    return data_base