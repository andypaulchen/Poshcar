# Compcar: Compare POSCAR, before and after relaxation

# Are you curious to find out how much atoms move before and after the relaxation step? This routine
# converts the direct coordinates to cartesian and then calculate the movement vector. I imagine you can
# plot the results by MATLAB and see which atoms move the furthest.

# Andy Paul Chen, Monday 20 April 2020, Cleveland Heights, Ohio

# Day #29 of lockdown. The archery target has arrived and I have been shooting some arrows. The punchbag
# loan from Krav Maga Warrensville has also come to great use. I go to Huaxin (North Randall) every
# Wednesday for groceries. The ladyboss laughs at me, and thinks I go there too often. "Have you eaten
# everything already?" she asks. Should I tell her I am only there to see her pretty face, and all that
# groceries are just a pretext? We had snow twice last week, which turned spring into winter briefly.
# I went out and saw Bambi eating fresh spring leaves by Doan Brook, but don't worry, I didn't shoot him.
# The big shale rock face at Frankie's Gate is always wet. I am wondering if this is what people call
# spring water.

from Cartesian import *

def data_atomsno(data1, data2):
    # Return list of elements and number of elements
    
    # Identify atomic species and number of atoms present in cell
    elem1 = re.findall(r'\w+', data1[atom_name_index].strip())
    num1 = re.findall(r'\d+', data1[atom_number_index].strip())
    elem2 = re.findall(r'\w+', data2[atom_name_index].strip())
    num2 = re.findall(r'\d+', data2[atom_number_index].strip())
    
    return (elem1==elem2) and (num1==num2)
    
def compcar(file1, file2, out_filename):
    # Read POSCAR data before (file1) and after (file2) relaxation, provide the movement vectors
    # under out_filename (cartesian coordinates). Output file is organized as columns of data,
    # as in POSCAR except there is a fourth column describing the distance moved by atom during relaxation
    # if Selective Dynamics is off in POSCAR, drift correction might need to be implemented
    # (not supported for now, probably not needed?)
    
    # Read input files
    data1 = readfile(file1)
    data2 = readfile(file2)
    
    if data_atomsno(data1, data2):
        print("Number of atoms by species: match!")
        # Execute main compcar routine
        dcindex = 8 if data_isSeldyn(data1) else 7 # Index of Direct/Cartesian line
        if not data_isCart(data1):
            data1 = data_switchCart(data1)
        if not data_isCart(data2):
            data2 = data_switchCart(data2)
        dataout = data2[:7]
        compflag = "Compcar\n"
        dataout += compflag
        
        for line in range(len(data1)):
            if line > dcindex:
                # Read
                comp1 = re.findall(r"-?\d+\.\d+", data1[line].strip())
                comp2 = re.findall(r"-?\d+\.\d+", data2[line].strip())
                # Compare
                co = [0,0,0]
                for i in range(3):
                    co[i] = float(comp2[i]) - float(comp1[i])
                # Fourth Column
                ab = math.sqrt(co[0]**2 + co[1]**2 + co[2]**2)
                # Write
                writeline = ls + "{:11f}".format(co[0]) + ls + "{:11f}".format(co[1]) + ls + "{:11f}".format(co[2]) + ls + "disp (angst.) =" + "{:11f}".format(ab) + "\n"
                dataout += writeline
                
        # write file
        writefile(dataout, out_filename)
    else:
        print("They are not the same system! No output file written.")
    