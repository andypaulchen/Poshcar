# atomsub: Algorithm to perform substitutional defects by doping (multiple atoms or a vacancy)

# (Atomsub 1 description):
# [This algorithm on the VASP POTCAR file creates a substitutional defect where a dopant atom occupied
# a native vacancy defect. This is done by removing a line from the POSCAR (vacancy) and appending said
# line at the end (substitutional doping). The list of elements (line 5) and number of elements per
# atom (line 6) are updated.]
# Update: This is a rewrite. I actually need to dope 2 substitutional atoms in my work. Stupid!
# The substituted atom index has been replaced by subatoms, a list of indices.

# Andy Paul Chen, Tuesday 31 March 2020, Cleveland Heights, Ohio

# Day #9 of lockdown. I called some friends back in Singapore, and showed them the new archery set that
# I bought for self-defense in the event of a societal collapse. They laughed at me for thinking to
# bring a bow and arrow to a gunfight. Gah! I'll show 'em.

# Andy Paul Chen, Friday 4 December 2020, Cleveland Heights, Ohio

# Hi, I'm back and I want to do vacancies now. It does not seem to be a functionality yet. What!
# Day ????? of Lockdown. I don't think the lockdown is a thing anymore. However some of my friends
# have da Covids. They are lucky not to be in hospital. As for myself I am 50% sure I do not have it.
# It seems no longer politically correct to say that the virus came from Wuhan. I am sad because it took
# this nice town where my mother went to school out of the spotlight unfairly. News articles post some
# update on lab results and then a few days later another article completely contradicting the first.
# This I attribute to the bad habit journalists have in using arXiv papers. I plan to spend Christmas
# reading a Dostoevsky, instead of flying home, as per custom.

# Andy Paul Chen, Monday 9 May 2022, Singapore

# I split up the function to a data_-headed thingimmabob. What's new in my life?
# I now have a Krav Maga practioner certificate, for what it's worth

from poshcar.seldyn import * # selective dynamics package

def atomsub(data, dopant_name, subatoms, verbose = True):
    # Print Input file information
    if verbose:
        print("List of elements (clean cell):"+ longspace + data[atom_name_index].strip())
        print("Number of atoms per element:" + longspace + data[atom_number_index].strip())
        print(horizont)
    
    # In the case of Selective Dynamics, the first coordinates are on line 8
    # Otherwise, it's 7
    addindex = 8 if is_seldyn(data) else 7
    
    # Read string of element symbols and number of atoms per element
    elem_list = re.findall(r'\w+', data[atom_name_index].strip())
    temp = re.findall(r'\d+', data[atom_number_index].strip())
    num_list = list(map(int, temp))
    
    # Remove substituted atoms from list (number-=len(subatoms))
    for k in range(len(subatoms)):
        modlo = subatoms[k]
        if modlo > sum(num_list):
            print("ERROR: One of the atoms does not exist!\n")
            break
        else:
            for j in range(len(num_list)):
                if modlo <= num_list[j]:
                    subd_element = j
                    num_list[subd_element] -= 1
                    s = subatoms[k] + addindex
                    switchout = data[s] # vacancy coordinates
                    del data[s]
                    if dopant_name in periodic_table:
                        data = data + [switchout]
                    break
                else:
                    modlo -= num_list[j]
                                  
    # Append symbol of dopant atom to list of elements
    # Append number of dopant atom to list of numbers
    if dopant_name not in periodic_table: dopant_name = ""
    data[atom_name_index] = ls + data[atom_name_index].strip() + ls + dopant_name + "\n"
    if verbose: print("List of elements (with dopant):", data[atom_name_index].strip())
    data[atom_number_index] = longspace
    for m in range(len(num_list)):
        data[atom_number_index] = data[atom_number_index] + str(num_list[m]) + ls
    if dopant_name not in periodic_table: dopant_number = ""
    else: dopant_number = str(len(subatoms))
    data[atom_number_index] = data[atom_number_index] + dopant_number + "\n"
    if verbose: print("Number of atoms per element: " + data[atom_number_index].strip())
    
    return data