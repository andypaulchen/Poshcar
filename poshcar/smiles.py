# smiles: Generate a 0D molecule cell of organic molecule from SMILES string
# Here we use the rdkit package as an aid. They got there first!!!!

# Andy Paul Chen, 25 November 2023

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from scipy.spatial.transform import Rotation as R
from poshcar.cartesian import *

def build_molecule(smiles):
    # From SMILES string, build a VASP cell containing an organic molecule
    
    my_mol = Chem.MolFromSmiles(smiles)
    my_mol_with_H=Chem.AddHs(my_mol)
    AllChem.EmbedMolecule(my_mol_with_H)
    AllChem.MMFFOptimizeMolecule(my_mol_with_H) # I suppose this is the MD line?
    MolBlock = Chem.MolToMolBlock(my_mol_with_H)  # contains coordinates
    MolBlockLines = MolBlock.splitlines()
    N = int(MolBlockLines[3].split()[0]) # total number of atoms
    
    # Read all data from the MolBlock!
    MBcoords = []
    MBelems = []
    for i in range(N):
        thiscoord = re.findall(r"-?\d+\.\d+", MolBlockLines[i+4].strip())
        thiscoord = [float(x) for x in thiscoord]
        MBcoords.append(thiscoord)
        MBelems.append([i, MolBlockLines[i+4].split()[3]])
            
    # Sort MBelems by element, putting C and H at head
    MBelems.sort(key=lambda x: x[1])
    allatomlist = [x[1] for x in MBelems]
    elemlist = list(set(allatomlist)) # enumerate unique atoms
    elemlist.sort()
    elemnos = [[u,allatomlist.count(u)] for u in elemlist]
    elemnos.sort(key=lambda x: x[0])
    elemnos = [str(x[1]) for x in elemnos]
    
    # Decide the size of box to put molecule in
    maxradius = max(-np.amin(MBcoords),np.amax(MBcoords))
    a = 2*maxradius + 10 # give margin of 10 angstroms between molecules, at least
    
    # Writing the lines
    astr = "{:11f}".format(a)
    data = [smiles+"\n", "1.0\n", 
            ls+astr+"  0.0  0.0\n", ls+"0.0  "+astr+"  0.0\n", ls+"0.0  0.0  "+astr+"\n",
            ls + ls.join(elemlist) + "\n", ls + ls.join(elemnos) + "\n", "Cartesian\n"]
    for j in MBelems:
        write = MBcoords[j[0]]
        data.append(ls + flpr.format(write[0]) + ls + flpr.format(write[1]) + ls + flpr.format(write[2]) + "\n")
    
    return data

def rotate_molecule(data, axis, angle):
    # Note: angle is in degrees, while axis is in 'x/y/z' or any combination thereof
    rotmatrix = R.from_euler(axis, angle, degrees=True)
    rotmatrix = rotmatrix.as_matrix()
    # See scipy documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.from_euler.html
    # Convert to Cartesian coordinates
    if not is_cart(data): data = switchcart(data)
    if is_seldyn(data): data = seldynswitch(data)
    dcindex = 8 if is_seldyn(data) else 7 # Index of Direct/Cartesian line
    for line in range(len(data)):
        if line > dcindex:
            coords = np.array(re.findall(r"-?\d+\.\d+", data[line].strip()))
            coords = coords.astype(float)
            newc = np.matmul(rotmatrix, coords)
            data[line] = ls + str(newc[0]) + ls + str(newc[1]) + ls + str(newc[2]) + "\n"
    return data

def center_molecule(data):
    B = basis(data) # Read lattice vectors
    vector = (B[0] + B[1] + B[2])/2
    data = translate(data, vector)
    return data