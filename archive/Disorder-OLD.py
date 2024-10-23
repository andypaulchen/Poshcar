# Disorder.py: Encoding site-disorder in a ".vasp" format
# Since the format ignores any term after the coordinate line (UNLESS selective dynamics is switched on),
# We can include atomic species and occupancy terms after them. The result is still readable.
# Note: do not use with selective dynamics!!!!

# 8 May 2024: I starts to implement the Atomic Sudoku code

from supercell import *
from atomsub import *
from Distance import *
from Additive import *
import random
import os

import ase.io.cif as asecif
from chgnet.model.model import CHGNet
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation

def desymmetrize(ciffile):
    # Reads in a .cif file, output another .cif file desymmetrised
    # Not used -- pymatgen doesn't handle this very well
    filename_header = ciffile[:-4] #output to same place as input
    # Step 1: Read the CIF file
    structure = Structure.from_file(ciffile)
    # Step 2: Symmetry reduction to P1 (primitive structure) while preserving partial occupancies
    sga = SpacegroupAnalyzer(structure, symprec=0.01)
    p1_structure = sga.get_primitive_standard_structure()
    # Step 3: Save the P1 structure with partial occupancies to a CIF file
    p1_structure.to(fmt="cif", filename=filename_header+'_P1.cif')


def cif2vasp_occ(ciffile, verbose = True):
    # Reads in a .cif file, output .vasp file with fractional occupancy information which mimics a Selective Dynamics file
    #desymmetrize(ciffile)
    filename_header = ciffile[:-4] #output to same place as input
    cifdata = asecif.read_cif(ciffile, index=-1, fractional_occupancies = True)
    
    #print(cifdata.__dict__['info']); cif data extractions
    occdict = cifdata.__dict__.get('info').get('occupancy') # dictionary containing species and occupancies
    sites = cifdata.__dict__.get('arrays').get('spacegroup_kinds') # array of sites, in sequence
    positions = cifdata.__dict__.get('arrays').get('positions') # array of site coordinates

    # write a vasp file here with ase, then read from it (yeah man that's dumb I know)
    cifdata.write(filename_header+'.vasp', format = 'vasp')
    vaspdata = readfile(filename_header+'.vasp')
    vaspdata.insert(7, "Site-Disordered Structure\n") # it starts with 'S' so VESTA can read this (so lame)
    vaspdata = vaspdata[:10] # keep header + 1 atom
    vaspdata[5] = " Ph \n"
    vaspdata[6] = "  1 \n" # Placeholder atom "Ph"

    cursor = 9
    for i in occdict:
        if int(i) in sites: j = int(i)
        else: j = max([x for x in sites if x < int(i)], default=None) 
        k = np.where(sites == j)[0][0] # extract coordinates
        toadd = next(iter(occdict[i].items())) # extract atoms
        vaspdata = addatom(vaspdata, toadd[0], positions[k], verbose = False)
        cursor += 1
        vaspdata[cursor] = vaspdata[cursor].rstrip() + ls + (toadd[0]+str(j)).ljust(5) + "  " + str(toadd[1]) + "\n"
    
    vaspdata = removeatom(vaspdata, [1], verbose = False) # Remove the placeholder atom
    # write data to vasp file again
    writefile(vaspdata, filename_header+'.vasp', verbose)
    return vaspdata

def composition(data, verbose = True):
    # Return vector containing fractional 
    df = elemindices(data, verbose)
    total = df['Occupancy'].sum()
    sum_by_spp = df.groupby('Species')['Occupancy'].sum()
    proportions = sum_by_spp/total
    if verbose: display(proportions)
    return proportions.values

def random_fill(data, crashmat, verbose = True):
    # Random fill of atomic sites
    # Input: disordered vasp poscar cell
    # Delete lines where occupancy is less than 1, and a dice-roll has decreed that it is removed
    # 24 Jul 2024 update: changed conformerID to string - too many permutations to store as int!!! (removed)
    # crashmat = crash matrix

    # Initializations
    backline = len(data)-1
    deletelist = [] # atoms to remove without cointoss

    while backline > 8:
        atomid = backline-8 # current atom index 
        occ = float(re.findall(r'\S+', data[backline].strip())[4]) # probability atom stays
        if occ < 1.0:
            if atomid in deletelist:
                data = atomsub(data,"vac",[atomid], verbose = False) # delete atom by crash
                if verbose: print("delete atom #", atomid, " by crash")
            elif random.random() > occ: 
                data = atomsub(data,"vac",[atomid], verbose = False) # delete atom by cointoss
                if verbose: print("delete atom #", atomid, "by cointoss")
                # need to increase occ for atoms that crash into them
            else:
                for x in range(1,atomid):
                    if crashmat[atomid-1][x-1] == 1: deletelist += [x]
        backline -= 1


def VirtualLibraryRandomFill(vaspfile, header, N, supercellsize, verbose = True, bond_threshold = 0.1, crash_threshold = 0.3, chgnet_on = False, targetfile = "."):
    # Source vaspfile is required to be a "disordered" vasp file
    # Generate a list of N virtual cells (text line-list format)
    # Check for target composition
    # supercellsize = array with 3 integers
    # header = folder name or path
    
    datalist = [] # init array
    complist = [] # composition
    cdlist   = [] # distance from target composition (Frobenius)
    bondlist = [] # bonding matrices (flattened)
    bdlist   = [] # distance from target bonding profile
    E_list   = [] # chgnet total energies

    data = readfile(vaspfile)
    try: os.mkdir('_operations/'+header)
    except: print('Folder already created')
    if chgnet_on: chgnet = CHGNet.load()
    
    if is_seldyn(data) and data[7].lower() == 'site-disordered structure\n': 
        atomoccs = elemindices(data, verbose)
        original_comp = composition(data, verbose) # register real compositional data
        partial_sites = (atomoccs['Occupancy'] < 1.0).sum()*np.prod(supercellsize)
        print("Number of partially-occupied sites: ", partial_sites)

        # Make target bonding matrix
        try: target = readfile(targetfile)
        except: target = data
        runiq, bme_targ, bma_targ = matrix_bonding_average(target, 'element', bond_threshold, verbose)

        # Making the supercell
        scz = supercellsize
        sc = supercell(data, scz[0], scz[1], scz[2]) # supercell
        i = 0 # initialize
        crashmat = crashtest(sc, crash_threshold, verbose = False)

        # Make virtual cells
        while i < N:
            print("i=",i)
            sci = sc.copy()
            sci = random_fill(sci, crashmat, verbose) # manipulate a copy of the supercell
            runiq, bme, bma = matrix_bonding_average(sci, 'element', bond_threshold, verbose)

            datalist.append(sci)
            #bondlist.append(bme.flatten())
            bondlist.append(bme)
            bdlist.append(np.linalg.norm(bme - bme_targ, 'fro')) # Bonding matrix differences / Frobenius Norm
            virtual_comp = composition(sci, verbose)
            complist.append(virtual_comp)
            try: compdist = distance(original_comp, virtual_comp)
            except: 
                print("Elements mismatch! (possibly all atoms of certain species removed?")
                compdist = "err" # this is not very probable just because of compounded probabilities!
            cdlist.append(compdist)
            writefile(sci, '_operations/'+header+'/'+str(i)+'.vasp')

            if chgnet_on:
                chgnet_structure = Structure.from_file('_operations/'+header+'/'+str(i)+'.vasp')
                prediction = chgnet.predict_structure(chgnet_structure)
                totalenergy = float(prediction['e'])
                if verbose: print(f"CHGNet-predicted energy (eV/atom)):\n{totalenergy}\n")
                E_list.append(totalenergy)
            i += 1

        if chgnet_on: summary = pd.DataFrame({'BondingProfile': bondlist, 'BondDiff' : bdlist, 'Composition': complist, 'CompDiff': cdlist, 'FormationEnergy': E_list})
        else: summary = pd.DataFrame({'BondingProfile': bondlist, 'BondDiff' : bdlist, 'Composition': complist, 'CompDist': cdlist})

        summary.to_csv('_operations/'+header+'/VirtualLibrary.csv')
        display(summary)
        print("Summary of results written to _operations/VirtualLibrary.csv")
        return datalist, summary
    else: print("This is not a site-disordered structure, dude!\n")
                    
def VirtualLibraryCif(ciffile, header, supercellsize, verbose = True, bond_threshold = 0.1, chgnet_on = False, targetfile = "."):
    # import disordered cell directly as cif
    # Check for target composition
    # supercellsize = array with 3 integers
    # header = folder name or path

    datalist = [] # init array
    complist = [] # composition
    cdlist   = [] # distance from target composition (Frobenius)
    bondlist = [] # bonding matrices (flattened)
    bdlist   = [] # distance from target bonding profile
    E_list   = [] # chgnet total energies

    # Prepare target subfolder and CHGNET
    try: os.mkdir('_operations/'+header)
    except: print('Folder already created')
    if chgnet_on: chgnet = CHGNet.load()

    # Make target bonding matrix
    vaspdata = cif2vasp_occ(ciffile, verbose = False)
    original_comp = composition(vaspdata, verbose)
    try: target = readfile(targetfile)
    except: target = vaspdata
    runiq, bme_targ, bma_targ = matrix_bonding_average(target, 'element', bond_threshold, verbose)

    # Load file into a structure, make supercell
    struct = Structure.from_file(ciffile)
    sup = struct.make_supercell(supercellsize, in_place = False)

    # Generate the cells
    to_virtual = OrderDisorderedStructureTransformation(no_oxi_states = True)
    resultat = to_virtual.apply_transformation(sup, return_ranked_list = N)
    print("Hear ye! size of dataset is ", len(resultat))

    # Write the structure to a VASP POSCAR file
    i = 0
    for entry in resultat:
        serial = str(i)
        vaspfilepath = '_operations/'+header+'/'+serial+'.vasp'
        entry['structure'].to(fmt="POSCAR", filename=vaspfilepath)
        entry_in_vasp = readfile(vaspfilepath, verbose)
        datalist.append(entry_in_vasp)

        # Bondmatrix block
        runiq, bme, bma = matrix_bonding_average(entry_in_vasp, 'element', bond_threshold, verbose)
        bondlist.append(bme)
        bdlist.append(np.linalg.norm(bme - bme_targ, 'fro')) # Bonding matrix differences / Frobenius Norm

        # Composition block
        virtual_comp = composition(entry_in_vasp, verbose)
        complist.append(virtual_comp)
        try: compdist = distance(original_comp, virtual_comp)
        except: 
            print("Elements mismatch! (possibly all atoms of certain species removed?")
            compdist = "err" # this is not very probable just because of compounded probabilities!
        cdlist.append(compdist)

        # CHGNET block
        if chgnet_on:
            chgnet_structure = Structure.from_file(vaspfilepath)
            prediction = chgnet.predict_structure(chgnet_structure)
            totalenergy = float(prediction['e'])
            if verbose: print(f"CHGNet-predicted energy (eV/atom)):\n{totalenergy}\n")
            E_list.append(totalenergy)
        i += 1

        if chgnet_on: summary = pd.DataFrame({'BondingProfile': bondlist, 'BondDiff' : bdlist, 'Composition': complist, 'CompDiff': cdlist, 'FormationEnergy': E_list})
        else: summary = pd.DataFrame({'BondingProfile': bondlist, 'BondDiff' : bdlist, 'Composition': complist, 'CompDist': cdlist})

    summary.to_csv('_operations/'+header+'/VirtualLibrary.csv')
    display(summary)
    print("Summary of results written to _operations/VirtualLibrary.csv")
    return datalist, summary

def VirtualLibrary(path, target, pauling_weight = 1, bond_threshold = 0.1, verbose = True):
    # Source: a folder (path) already populated by virtual cells from supercell software
    # No further filling required!!
    # Check for target composition (target is compulsory)
    # header = folder name or path
    
    datalist = [] # init array
    complist = [] # composition
    cdlist   = [] # distance from target composition (Frobenius)
    bondlist = [] # bonding matrices
    bdlist   = [] # distance from target bonding profile
    E_list   = [] # chgnet total energies

    # Load CHGNET bc we defo need to use this
    chgnet = CHGNet.load()

    # Take care of the target file (could be source file input to supercell)
    # or maybe even a user-defined cell for correlated disorder
    runiq, bme_targ, unav = matrix_bonding_average(target, 'element', bond_threshold, verbose = True)
    original_comp = composition(target, verbose) # register real compositional data

    # Iterate over all .cif files in folder
    for ciffile in os.listdir(path):
        if ciffile.endswith('.cif'):
            filename_header = ciffile[:-4] #output to same place as input

            # Convert cif to vasp for our matrices
            cifdata = Structure.from_file(path+ciffile)
            vaspfilepath = path + filename_header + '.vasp'
            cifdata.to(filename = vaspfilepath, fmt = "poscar")
            virtual = readfile(vaspfilepath)
            if verbose: printvaspdata(virtual)
            runiq, bme, unaverageness = matrix_bonding_average(virtual, 'element', bond_threshold, bme_correlated = bme_targ, verbose = verbose)

            # append to datalist
            datalist.append(virtual)

            # append to bond lists
            bondlist.append(bme)
            # bonding matrix distance = difference in bme + Pauling-5 penalty
            bd_metric = np.linalg.norm(bme - bme_targ, 'fro') 
            print("Bonding matrix distance: ", bd_metric, " Unaverageness: ", unaverageness)
            bdlist.append(bd_metric + pauling_weight*unaverageness)
            #bdlist.append(bd_metric)

            # append to composition lists
            virtual_comp = composition(virtual, verbose)
            complist.append(virtual_comp)
            try: compdist = distance(original_comp, virtual_comp)
            except: 
                print("Elements mismatch! (possibly all atoms of certain species removed?")
                compdist = "err" # this is not very probable just because of compounded probabilities!
            cdlist.append(compdist)

            # perform chgnet and append total energy to list
            chgnet_structure = Structure.from_file(vaspfilepath)
            prediction = chgnet.predict_structure(chgnet_structure)
            totalenergy = float(prediction['e'])
            if verbose: print(f"CHGNet-predicted energy (eV/atom)):\n{totalenergy}\n")
            E_list.append(totalenergy)

    summary = pd.DataFrame({'BondingProfile': bondlist, 'BondDiff' : bdlist, 'Composition': complist, 'CompDiff': cdlist, 'FormationEnergy': E_list})
    summary.to_csv(path+'/VirtualLibrary.csv')
    display(summary)
    print("Summary of results written to _operations/VirtualLibrary.csv")
    return datalist, summary