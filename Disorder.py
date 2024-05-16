# Disorder.py: Encoding site-disorder in a ".vasp" format
# Since the format ignores any term after the coordinate line (UNLESS selective dynamics is switched on),
# We can include atomic species and occupancy terms after them. The result is still readable.
# Note: do not use with selective dynamics!!!!

# 8 May 2024: I starts to implement the Atomic Sudoku code

from Supercell import *
from AtomSub import *
from Distance import *
import ase.io.cif as asecif
import ase.io.vasp as asevasp
from ase.geometry import get_duplicate_atoms
import random
import os

def cif2vasp_occ(ciffile):
    # Reads in a .cif file, output .vasp file with fractional occupancy information which mimics a Selective Dynamics file
    filename_header = ciffile[:-4] #output to same place as input
    cifdata = asecif.read_cif(ciffile, fractional_occupancies = True)

    #print(cifdata.__dict__['info'])
    occdict = cifdata.__dict__.get('info').get('occupancy') # dictionary containing species and occupancies
    sites = cifdata.__dict__.get('arrays').get('spacegroup_kinds') # array of sites, in sequence

    # write a vasp file here with ase, then read from it (yeah man that's dumb I know)
    cifdata.write(filename_header+'.vasp', format = 'vasp')
    vaspdata = readfile(filename_header+'.vasp')
    vaspdata.insert(7, "Site-Disordered Structure\n") # it starts with 'S' so VESTA can read this (so lame)

    for line in range(len(vaspdata)):
        if line > 8: 
            ordinal = line-9
            site_id = str(sites[ordinal]) # Which site the atom is on
            site_spp_occ = occdict.get(site_id) # key is species, item is occupancy
            spp = list(site_spp_occ.keys())[0]
            occ = site_spp_occ.get(spp)
            vaspdata[line] = vaspdata[line].rstrip() + ls + (spp+site_id).ljust(5) + "  " + str(occ) + "\n"

    # write data to vasp file again
    writefile(vaspdata, filename_header+'.vasp')

def composition(data):
    # Return vector containing fractional 
    df = ElemIndices(data)
    total = df['Occupancy'].sum()
    sum_by_spp = df.groupby('Species')['Occupancy'].sum()
    proportions = sum_by_spp/total
    display(proportions)
    return proportions.values

def disordered_to_virtual(data):
    # Random fill of atomic sites
    # Delete lines where occupancy is less than 1, and a dice-roll has decreed that it is removed
    
    crash_threshold = 0.4 # Set the threshold here!
    # Initializations
    conformerID = 0 # Constructing the ConformerID with the binary method
    cid_place = 0
    backline = len(data)-1
    
    while backline > 8: 
        occ = float(re.findall(r'\S+', data[backline].strip())[4]) # probability atom stays
        if occ < 1.0:
            if random.random() > occ: data = AtomSub(data,"vac",[backline-8]) # delete atom
            else: conformerID += 2**cid_place
            cid_place += 1
        backline -= 1
        
    # Check for crashing atoms, target composition
    if Crashtest(data, crash_threshold): 
        print("Cell accepted.")
        crashtest_passed = True
        data = SeldynSwitch(data) # Switch for virtual cell
    else: 
        print("Rejected: Atoms too close together.")
        crashtest_passed = False
    return data, conformerID, crashtest_passed

def VirtualLibrary(vaspfile, header, N, supercellsize):
    # Source vaspfile is required to be a "disordered" vasp file
    # Generate a list of N virtual cells (text line-list format)
    # Check for target composition
    # Supercellsize = array with 3 integers
    # header = folder name
    
    datalist = [] # init array
    cidlist = [] # conformerID list
    crashlist = [] # pass/fail grade for crashtest
    complist = [] # composition-distance
    data = readfile(vaspfile)
    try: os.mkdir('_operations/'+header)
    except: print('Folder already created')
    
    if isSeldyn(data) and data[7].lower() == 'site-disordered structure\n': 
        # Set up conformerID and related conditions related to combinatorics
        atomoccs = ElemIndices(data)
        original_comp = composition(data) # register real compositional data
        partial_sites = (atomoccs['Occupancy'] < 1.0).sum()
        print("Number of partially-occupied sites: ", partial_sites)
        total_possible_conform = 2**partial_sites
        print("Number of possible conformers: ", total_possible_conform)
        N1 = min(N, total_possible_conform) # limit to the smaller number of conformers
    
        scz = supercellsize
        sc = Supercell(data, scz[0], scz[1], scz[2]) # supercell
        i = 0 # initialize
        while i < N1:
            print("i=",i)
            sci = sc.copy()
            sci, cid, crash = disordered_to_virtual(sci) # manipulate a copy of the supercell
            if cid not in cidlist:
                datalist.append(sci)
                cidlist.append(cid)
                crashlist.append(crash)
                virtual_comp = composition(sci)
                try: compdist = distance(original_comp, virtual_comp)
                except: 
                    print("Elements mismatch! (possibly all atoms of certain species removed?")
                    compdist = "err"
                complist.append(compdist)
                writefile(sci, '_operations/'+header+'/'+str(cid)+'.vasp')
                i = i + 1
        summary = pd.DataFrame({'ConformerID': cidlist, 'Passed Crashtest': crashlist, 'CompDist': complist})
        summary.to_csv('_operations/'+header+'/VirtualLibrary.csv')
        display(summary)
        print("Summary of results written to _operations/VirtualLibrary.csv")
        return datalist
    else: print("This is not a site-disordered structure, dude!\n")
                    

    