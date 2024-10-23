# cquery: Easy look-up of atom coordinates

# Listen. You put in number x. Program find x-th atom on list, give coordinate. Давайте?
# Yuo can do the same in VESTA. Unfortunately, clicking in VESTA can be annoying. 
# This can be used as a double-checking apparatus before using other functions,
# in cases where look-up by index is faster (it does happen sometimes)

# Andy Paul Chen, 26 March 2020 (Coronavirus Season!), Cleveland Heights, Ohio

# Andy Paul Chen, Monday, 9 August 2021, Little Italy, Cleveland, Ohio (National Day of Singapore)
# data_cQuery dadded

# Changed name from Coords to cquery on 7 November 2023

from poshcar.seldyn import * # selective dynamics package

def cquery(data, atomno, verbose = True):
    # Identify atomic species and number of atoms present in cell
    # Returns elem
    elem_list = re.findall(r'\w+', data[atom_name_index].strip())
    num_list = re.findall(r'\d+', data[atom_number_index].strip())
    
    # Index of atoms (translate from ordinal -> location in file)
    addindex = 8 if is_seldyn(data) else 7
    
    # Now to identify species and ordinal number!
    #Read string of numbers
    no_list = list(map(int, num_list))
    
    modlo = atomno
    if modlo > sum(no_list):
        print("ERROR: Atom does not exist!\n")
    else:
        for j in range(len(no_list)):
            if modlo <= no_list[j]:
                # ID element
                elem = elem_list[j]
                if verbose: print("Atom ID: ", atomno, " / ", elem, modlo)
                # Print coordinates
                if verbose: print("Coordinates: ", data[atomno+addindex])
                coords = data[atomno+addindex].split()
                coords = [float(part) for part in coords]
                break
            else:
                modlo -= no_list[j]
                
    return elem, coords