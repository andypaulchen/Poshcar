# Perovskite: Generates a perovskite POSCAR 

# Andy Paul Chen, Thursday 17 February 2022, Singapore

from Cartesian import * # SHMOSCAR core package
ls = longspace
lls = ls + "   "
Organics = ["FA", "MA"]

def data_perovABX3(a, A, B, X):
    data = [A+B+X+"3\n", "1.0\n",
            ls+"{:11f}".format(a)+"  0.0  0.0\n", ls+"0.0   "+"{:11f}".format(a)+"   0.0\n", ls+"0.0  0.0  "+"{:11f}".format(a)+"\n", 
            ls+X+ls+B+ls+A+"\n", ls+"3    1    1\n","Direct\n",
            ls+"0.5  0.0  0.5\n",ls+"0.5  0.5  0.0\n",ls+"0.0  0.5  0.5\n",ls+"0.5  0.5  0.5\n",ls+"0.0  0.0  0.0\n"]
    if (A == "MA"):
        data[0] = "[MA](CH3-NH3)"+B+X+"3\n"
        data[5] = ls+X+ls+B+ls+"C"+ls+"N"+ls+"H\n"
        data[6] = ls+"3    1    1    1    6\n"
        data = data_switchCart(data)
        data[12] = lls+"0.585428"+lls+"-0.02016"+lls+"0.084438\n"
        data.extend([lls+"-0.88338"+lls+"-0.02016"+lls+"-0.18486\n",
                     lls+"-1.35278"+lls+"0.815464"+lls+"0.209024\n", 
                     lls+"-1.10662"+lls+"-0.02016"+lls+"-1.20319\n", 
                     lls+"1.022068"+lls+"0.877012"+lls+"-0.36352\n", 
                     lls+"0.752968"+lls+"-0.02016"+lls+"1.167930\n", 
                     lls+"1.022068"+lls+"-0.91733"+lls+"-0.36352\n", 
                     lls+"-1.35278"+lls+"-0.85578"+lls+"0.209024\n"])
    if (A == "FA"):
        data[0] = "[FA](H2N-CH-NH2)"+B+X+"3\n"
        data[5] = ls+X+ls+B+ls+"C"+ls+"N"+ls+"H\n"
        data[6] = ls+"3    1    1    2    5\n"
        data = data_switchCart(data)
        data.extend([lls+"-0.59155"+lls+"0.000000"+lls+"-1.17165\n",
                     lls+"-0.59155"+lls+"0.000000"+lls+"1.171655\n",
                     lls+"-1.60606"+lls+"0.000000"+lls+"1.309395\n",
                     lls+"-1.60606"+lls+"0.000000"+lls+"-1.30940\n",
                     lls+"1.096248"+lls+"0.000000"+lls+"0.000000\n",
                     lls+"-0.00878"+lls+"0.000000"+lls+"-2.00780\n",
                     lls+"-0.00878"+lls+"0.000000"+lls+"2.007801\n"])
    return data

def perovABX3(out_filename, a, A, B, X):
    # given lattice parameter a (in Angstroms), construct cell of ABX3 and write in out_filename
    # Order of atoms: X, B, A (A can be a set of atoms)
    
    if ((A in Periodic_Table) or (A in Organics)) and (B in Periodic_Table) and (X in Periodic_Table):
        dataout = data_perovABX3(a, A, B, X)
        writefile(dataout, out_filename)
    else:
        print("Error: Unidentifiable elements!\n")
    