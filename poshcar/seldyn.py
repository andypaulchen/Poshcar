# Seldyn.py: package for manipulating selective dynamics parameters in POSCAR

# Andy Paul Chen, 28 March 2020 (Lockdown Day #5)

# I have friends who might have contracted Covid. They have self-isolated. One of them, a very
# good friend of mine, lost his sense of smell and taste, but later got it back.
# May it stay that way!

# I am well save for a weird feeling in my soft palate. My body seems to be trying (but not
# very hard) to have a cold.
# I bought a bow and forearm guard for defense in the event of a societal meltdown. The arrows
# will come soon, hopefully.

from poshcar import *

def is_seldyn(data):
    # Read list of lines extracted from file (data) and determine if selective dynamics
    # is switched on or not.
    # True if yes, false if no
    # Sehr einfach!
    
    # Index of atoms (translate from ordinal -> location in file)
    if data[7][0].upper() == 'S': return True
    else: return False

def seldynswitch(data):
    # if data has no selective dynamics, turn it on
    # and switch all flags to false (F). If selective dynamics is on, remove the line containing
    # 'selective dynamics' and remove all flags.

    # Whittle the tail (CONTCAR)
    dcindex = 8 if is_seldyn(data) else 7 # Index of Direct/Cartesian line
    num = re.findall(r'\d+', data[atom_number_index].strip())
    snum = sum(list(map(int,num)))
    data = data[:(snum + dcindex + 1)] # whittle the tail
    
    # Test is_seldyn condition, remove line if true, add line if false
    # Remove flags if true, add flags (F F F) if false
    if is_seldyn(data):
        for line in range(len(data)):
            if line > 8:
                coords = re.findall(r"-?\d+\.\d+", data[line].strip())
                data[line] = ls + str(coords[0]) + ls + str(coords[1]) + ls + str(coords[2]) + "\n"
        del data[7]
    else:
        data.insert(7, "Selective Dynamics\n")
        for line in range(len(data)):
            if line > 8: data[line] = data[line].rstrip() + ls + "F F F\n"
    return data
    
# Andy Paul Chen, 1 April 2020 (Lockdown Day #9)

# Hi, I'm back. My soft palate is good now, but there is a stupid cough. I want to write
# some more functions. The "quiver" for the arrows has been delivered to my place. It is
# exactly the same thing as that tube you use to carry conference posters around.
# 5.40pm: Now the arrows are also here

def setflags(data, TF, setatoms):
    # TF: any combination of "(T/F)(T/F)(T/F)" as flags
    # e.g. set TF = 'F F T' to relax atom in +z direction only
    # setatoms: [list] of atoms affected by operation
    
    # Case-insensitive rigging
    TF = TF.upper()
    LTF = list(TF)
    TF_valid = len(TF)==3 and LTF.count('T')+LTF.count('F')==3
    if not TF_valid: print("Error in argument 3: 3 instances of 'T' or 'F' expected!")
        
    # If selective dynamics not switched on, do so and set all flags to F
    if not is_seldyn(data):
        print("Selective dynamics not switched on in POSCAR. Well dang, I'll switch it anyway...")
        data = seldynswitch(data)
        print("All flags set to F")
    
    # Read string of element symbols and number of atoms per element
    atoms_valid = True
    temp = re.findall(r'\d+', data[atom_number_index].strip())
    num_list = list(map(int, temp))
    if max(setatoms) > sum(num_list):
        print("Error in argument 4: invalid atom index!")
        atoms_valid = False
    else:
        for index in setatoms:
            index += 8
            data[index] = (data[index].rstrip())[:-6]
            data[index] = data[index].rstrip() + longspace + LTF[0] + " " + LTF[1] + " " + LTF[2] + "\n"
    
    return data