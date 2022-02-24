# Coords: Easy look-up of atom coordinates

# Listen. You put in number x. Program find x-th atom on list, give coordinate. Давайте?
# Yuo can do the same in VESTA. Unfortunately, clicking in VESTA can be annoying. 
# This can be used as a double-checking apparatus before using other functions,
# in cases where look-up by index is faster (it does happen sometimes)

# Andy Paul Chen, 26 March 2020 (Coronavirus Season!), Cleveland Heights, Ohio

# Andy Paul Chen, Monday, 9 August 2021, Little Italy, Cleveland, Ohio (National Day of Singapore)
# data_cQuery dadded

from Seldyn import * # selective dynamics package

def data_cQuery(data,atomno):
    # Identify atomic species and number of atoms present in cell
    elem_list = re.findall(r'\w+', data[atom_name_index].strip())
    num_list = re.findall(r'\d+', data[atom_number_index].strip())
    
    # Index of atoms (translate from ordinal -> location in file)
    addindex = 8 if data_isSeldyn(data) else 7
    
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
                print("Atom ID: ", atomno, " / ", elem, modlo)
                # Print coordinates
                print("Coordinates: ", data[atomno+addindex])
                break
            else:
                modlo -= no_list[j]
    
def cQuery(in_filename, atomno):
    # Take (in_filename) file, prints coordinates line of (atomno)-th atom as string
    # Also prints element, ordinal number in list of atoms in element sublist (e.g. Al16)
    # and flags for selective dynamics, as an extra bonus!
    
    # Read input file
    data = readfile(in_filename)
    
    # Run data_cQuery
    data_cQuery(data,atomno)