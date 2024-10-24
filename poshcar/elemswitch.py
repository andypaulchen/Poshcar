# ElemSwitch: Algorithm to switch two elements in the list of elements in POSCAR

# This switches the order of elements in the POSCAR. Sometimes the visualisation software
# does it on its own without warning you. Isn't it annoying? Looking at you, VESTA

# Andy Paul Chen, 26 March 2020 (Coronavirus Season!), Cleveland Heights, Ohio

# Day #3 of lockdown. I had symptoms two days ago. I don't have them anymore, and am no
# longer fearing for my life. Jsyk, I'm not been playing with this for long enough to be
# absolutely sick of Python, so y'alls can infer that I'm green. Wish me luck.
# I do my documentation the best I can in the comments. You are very welcome!

from poshcar.seldyn import * # selective dynamics package

def elemswitch(data):
    # Switch elements based on input prompt
    
    # Print Input file information
    print("List of elements:"+ longspace + data[atom_name_index].strip())
    print("Number of Atoms:" + longspace + data[atom_number_index].strip())
    print(horizont)
    
    # Identify atomic species
    elem_list = re.findall(r'\w+', data[atom_name_index].strip())
    num_list = re.findall(r'\d+', data[atom_number_index].strip())

    # Index of atoms (translate from ordinal -> location in file)
    # In the case of Selective Dynamics, the first coordinates are on line 8
    # Otherwise, it's 7
    addindex = 8 if is_seldyn(data) else 7
    
    # Prompt user to say which ones to switch
    user = input("Which two elements do you want to switch?: ")
    user_list = re.findall(r'\w+', user.strip())
    print("Chosen atomic species to switch: ", user_list)
    
    if all(x in elem_list for x in user_list) and len(user_list) == 2:
        card = len(elem_list) # cardinality of set
        
        #subset test, require 2 elements
        order = [z for z in range(card)]
        for i in range(card):
            if (user_list[0] == elem_list[i]): sw1 = i
            if (user_list[1] == elem_list[i]): sw2 = i
        switchb = order[sw1]
        order[sw1] = order[sw2]
        order[sw2] = switchb
        
        # Switch elements in element name and number lines
        data[atom_name_index] = ""
        data[atom_number_index] = ""
        for k in range(card):
            data[atom_name_index] += (longspace + elem_list[order[k]])
            data[atom_number_index] += (longspace + str(num_list[order[k]]))
        data[atom_name_index] += "\n"
        data[atom_number_index] += "\n"
        print("New Element Order:", longspace, data[atom_name_index].strip())
        print("Number of Atoms:", longspace, data[atom_number_index].strip())
        
        # identify blocks (line # in file of last atom of each element)
        blox = [0]*(card+1)
        startl = 0
        endl = addindex + 1
        dataout = data[startl:endl] # includes only the header at this point
        for j in range(card):
            startl = endl
            endl = int(endl) + int(num_list[j])
            blox[j] = data[startl:endl]
    
        # Concatenate poscar file 
        for m in range(card): dataout += blox[order[m]]       
        return dataout
    else:
        print("ERROR: Invalid switching duplet! Output file not written")


# Day 37 of Lockdown: 28 April 2020. Time passes quickly like the arrows I now shoot gleefully
# at the target in the basement. I meet my goals also like the arrows meeting the target at the 
# basement, in other words badly.
        
def elemset(data, pos, Sp2, verbose = True):
    # Take in_filename, change all atom of species in ordinal position pos (1,2,3...etc) to Sp2
    # Print Input file information
    if verbose: print("List of elements:"+ longspace + data[atom_name_index].strip())
    
    # Identify atomic species
    elem_list = re.findall(r'\w+', data[atom_name_index].strip())
    print(len(elem_list))
    
    if type(pos)==int and pos > 0 and pos <= len(elem_list):
        if Sp2 in periodic_table:
            elem_list[pos-1] = Sp2
            # Construct new line
            data[atom_name_index] = ""
            for k in range(len(elem_list)):
                data[atom_name_index] += (longspace + elem_list[k])
            data[atom_name_index] += "\n"
            if verbose: print("New Element Order:", longspace, data[atom_name_index].strip())
            return data
        else: print("ERROR: Invalid element!")
    else: print("ERROR: Invalid parameter #2!")