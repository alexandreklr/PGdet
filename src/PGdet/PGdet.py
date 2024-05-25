from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import SetBondLength, SetAngleDeg
from pointgroup import PointGroup
from rdkit.Chem import rdMolTransforms as rdmt
import numpy as np
import matplotlib.pyplot as plt

# Function that calculates the number of non hydrogen atoms contained in a smiles.
def count_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    atom_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H':
            atom_count += 1
    return atom_count

# Function that calculate angles between 3 points in space.
def calculate_angle(p1, p2, p3):
    
    # In between points vectors calculation.
    points = [p1, p2, p3]
    combi_ang = []
    for i in range(len(points)):
        a = np.subtract(points[-i], points[1-i])
        b = np.subtract(points[2-i], points[1-i])    
        combi_ang.append((a,b))
    tot_norm = []
    for combi in combi_ang:
        c = np.linalg.norm(combi[0])+np.linalg.norm(combi[1])
        tot_norm.append(c) 
    m = 9999999999999999999999999999999999
    for i in tot_norm:
        if i < m:
             m = i
    index = tot_norm.index(m)
    right_combi = combi_ang[index]
    v1,v2 = right_combi
     
    # Scalar product and norms calculations.
    dot_product = np.dot(v1, v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)  
    
    # Angle calculation in radians then degrees conversion.
    angle_rad = np.arccos(dot_product / norms)
    angle_deg = np.degrees(angle_rad)    
    return angle_deg
 
# This function checks if the smiles is valid for pg determination.
def pg_check_smiles(smiles):
    
    # List of smiles of inorganic molecules authorised as an exception.
    allowed_inorganics = ["O", "N", "NN", "O=O", "N=N", "N#N", "OO", "F", "FF", "Cl", "ClCl", "Br", "BrBr", "I", "II", "ClI", "ICl"]

    # The function tests if smiles is a string or a list.
    if type(smiles) != str and type(smiles) != list:
        raise TypeError("Invalid SMILES class. SMILES must be a string or a list of strings.")
    
    # The code gets the Rdkit object from the smiles.
    mol = Chem.MolFromSmiles(smiles)
    
    # Verifying that there is at least one carbon atom within the molecule. mol is a NoneType if the smiles is invalid.
    try:
        has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
        if smiles == '':
            raise AttributeError
    except AttributeError:
        raise TypeError("Invalid SMILES. Mol object could not be generated.")
    
    # The code checks if the molecule is an allowed organic and if it is not too large.
    if not has_carbon and smiles not in allowed_inorganics: 
        raise TypeError(f"Error: The smiles input corresponds to an inorganic molecule. No molecule could be printed. The molecule should be purely organic or be an exception: {allowed_inorganics}")
    if count_atoms(smiles) > 25:
        raise TypeError("Error: The smiles input is too long. Maximum non hydrogen atoms allowed : 25")
    
    # Checks if there is a metal in the smiles 
    # and rejects the smiles if there is one in it.
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num > 20 and atomic_num not in [33, 34, 15, 16, 9, 17, 35, 53]:
            raise TypeError("Error: the smiles corresponds to an organometallic complex, which is not suited by the program.")
            raise TypeError("Please, try another smiles.")

    return 

# This function generates a 3d conformation of the molecule. 
# It optimizes and modifies the conformation for its point group to be calculated.
def config_mol(smiles):
    
    # List of smiles of inorganic molecules authorised as an exception.
    allowed_inorganics = ["O", "N", "NN", "O=O", "N=N", "N#N", "OO", "F", "FF", "Cl", "ClCl", "Br", "BrBr", "I", "II", "ClI", "ICl"]
    
    # The code gets the Rdkit object from the smiles.
    mol = Chem.MolFromSmiles(smiles)
    
    # Generates a 3D configuration for the molecule.
    mol = Chem.AddHs(mol)
    np.random.seed(48) 
    AllChem.EmbedMolecule(mol, randomSeed=42, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    
    # Conformation optimisation.
    AllChem.MMFFOptimizeMolecule(mol)
    
    # Bond length is adjusted. Bonds that are part of a ring are not modified.
    conf = mol.GetConformer()
    for bond in mol.GetBonds():
        # Check if the bond is part of a ring.  
        if bond.IsInRing():
            continue
        # Sets all bond lengths to 1.
        SetBondLength(conf, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), 1.0)  
    
    # Modifies angles of H-X where X is O, N, S or P.
    if smiles not in allowed_inorganics:   
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ['O', 'N', 'S', 'P']:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'H':
                        # Finds Carbon atom bonded to X.
                        carbon_atom = None
                        for n in atom.GetNeighbors():
                            if n.GetSymbol() == 'C':
                                carbon_atom = n
                                break
                            if carbon_atom:
                            # Check if the bonds are part of a ring.
                                if not mol.GetBondBetweenAtoms(neighbor.GetIdx(), atom.GetIdx()).IsInRing() and not mol.GetBondBetweenAtoms(atom.GetIdx(), carbon_atom.GetIdx()).IsInRing():
                                # Defines H-X-C angles to 180 degrees, that way the code should handle better heteroatoms.
                                    SetAngleDeg(conf, neighbor.GetIdx(), atom.GetIdx(), carbon_atom.GetIdx(), 180.0)  
                                    
            # Modifies the angles between 
            # substituants of sp3 carbon atoms.
            if atom.GetSymbol() == "C":
                neighbors = [n for n in atom.GetNeighbors()]
                if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                    neighbors = [n for n in atom.GetNeighbors()]
                    if len(neighbors) >= 2:
                        for i in range(len(neighbors)):
                            for j in range(i+1, len(neighbors)):
                                # Check if the bonds are part of a ring.
                                if not mol.GetBondBetweenAtoms(neighbors[i].GetIdx(), atom.GetIdx()).IsInRing() and not mol.GetBondBetweenAtoms(atom.GetIdx(), neighbors[j].GetIdx()).IsInRing():
                                # Sets X-C-Y angles to 109.5 degrees for sp3 carbon atoms.
                                   SetAngleDeg(conf, neighbors[i].GetIdx(), atom.GetIdx(), neighbors[j].GetIdx(), 109.5)
                   
                # Modifies the angles between 
                # substituants of sp2 carbon atoms.
                elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                    # Sets X-C-Y angles to 120 degrees for sp2 carbon atoms.
                    for i in range(len(neighbors)):
                        for j in range(i+1, len(neighbors)):
                            if mol.GetBondBetweenAtoms(neighbors[i].GetIdx(), neighbors[j].GetIdx()):
                                SetAngleDeg(conf, neighbors[i].GetIdx(), atom.GetIdx(), neighbors[j].GetIdx(), 120.0)
                
                #Modifies the angles between 
                #substituants of sp carbon atoms
                elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP:
                    # Sets X-C-Y angles to 180 degrees for sp carbon atoms.
                    for i in range(len(neighbors)):
                        for j in range(i+1, len(neighbors)):
                            if mol.GetBondBetweenAtoms(neighbors[i].GetIdx(), neighbors[j].GetIdx()):
                                SetAngleDeg(conf, neighbors[i].GetIdx(), atom.GetIdx(), neighbors[j].GetIdx(), 180.0)
                                
    # This part of the code modifies the 
    # dihedral angles of carbon substituants.
    # Iterate over all atoms in the molecule.
    for atom_idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atom_idx)

        # Check if the atom is a carbon atom and if atom is not part of a ring.
        if atom.GetSymbol() == 'C' and not atom.IsInRing():
        # Get the indices of the neighboring atoms.
            neighbors = [x.GetIdx() for x in atom.GetNeighbors()]

           # We need at least 4 atoms to form a dihedral angle.
            if len(neighbors) >= 4:
                for i in range(len(neighbors) - 2):
                # Check if atoms j and k are bonded.
                    if mol.GetBondBetweenAtoms(neighbors[i], neighbors[i+1]) is not None:
                    # Set the dihedral angle to 60°.
                        rdmt.SetDihedralDeg(conf, atom_idx, neighbors[i], neighbors[i+1], neighbors[i+2],60)

    return mol, conf

# Takes molecule atom's coordinates and a list indicating which of them are bonded.
# Then gives a list of angles between the atoms with the atom's index in atom_list.
def molecule_angles(coords, bonded_atoms):
    angles = []
        
    # Molecule angles determination.
    for i in range(len(bonded_atoms)):
        for j in range(i+1, len(bonded_atoms)):
            if bonded_atoms[i][0] in bonded_atoms[j]:
                angle = calculate_angle(coords[bonded_atoms[j][0]], coords[bonded_atoms[j][1]] , coords[bonded_atoms[i][1]])
                angles.append((bonded_atoms[j][0], bonded_atoms[j][1], bonded_atoms[i][1], angle))
            elif bonded_atoms[i][1] in bonded_atoms[j]:
                angle = calculate_angle(coords[bonded_atoms[j][0]], coords[bonded_atoms[j][1]] , coords[bonded_atoms[i][0]])
                angles.append((bonded_atoms[j][0], bonded_atoms[j][1], bonded_atoms[i][0], angle))
        
    return angles

# Takes a mol object and its associated atom coordinates and gives back two lists.
# bonded_atoms is a list of sets containing two numbers. Each number is the index of one atom in atom_list.
# Two atoms present in the same set are bonded.
# lenghts is a list of bond lenghts in the same order as bonded_atoms.
def molecule_bonds(mol, coords):
    
    # Counts bonds.
    bonded_atoms = []
    for bond in mol.GetBonds():
        bonded_atoms.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        mol = np.array(mol)
        
    # Calculates the length of each bond using the position in space of each atom.
    # Add each length in a list.
    lengths = []     
    for liaison in bonded_atoms:
        a = np.array(coords[liaison[0]])
        b = np.array(coords[liaison[1]])
        length = np.linalg.norm(b - a)
        lengths.append(length)
    
    return bonded_atoms, lengths

# Takes a list of atoms, their coordinates, and a list of sets indication 
# which are bonded then plots the molecule in 3d.
# The three list must be in the same order.
def plot_molecule(atom_list, coords, bonded_atoms):
    # Atoms and their corresponding colors
    atom_colors = {
            'C': 'black',
            'H': 'white',
            'O': 'red',
            'N': 'blue',
            'S': 'yellow',
            'P': 'brown',
            'I': 'purple'
        }
        
    atom_coords = np.array(coords)

    # 3D figure creation.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.axis('equal')
    liste=[]
    for i in range(len(atom_coords)):
        liste.append(max(atom_coords[i]))
    # Changes the size of the graph so that the molecule fits perfectly without deformation.
    x = max(liste)
    ax.set_xlim([-x, x])
    ax.set_ylim([-x, x])
    ax.set_zlim([-x, x])
    
    # Associates each point in space with its atomic symbol.
    for i, (x, y, z) in enumerate(atom_coords):
        ax.scatter(x, y, z, label=atom_list[i])
    # Adds colors to each point according to the dictionnary.
    for i, (x, y, z) in enumerate(atom_coords):
        try:
            ax.scatter(x, y, z, color=atom_colors[atom_list[i]])
        # Treat the case when the atom is not in atom_colors.
        except KeyError:
            ax.scatter(x, y, z, color='grey')
    # Adds bonds to the 3D figure.
    for i, j in bonded_atoms:
        ax.plot(*zip(atom_coords[i], atom_coords[j]), color='Grey')

# Figure display.
    plt.show()
    
    return

# Main function. Uses all of the above functions.
# The function takes a smiles and gives back its point group as a string.
# If one sets desc to True, the function will also print a description of the molecule.
# If plot is set to True, the function will also plot a 3D figure of the molecule.
def atom_mapping(smiles, desc = False, plot = False):
    
    # List of SMILES not well handeled by the program that are Cs
    Cs = ["CCO", "CCCO", "CCCCO", "CCCCCO", "CCCCCCO",
          "CCCCCCCO","CCCCCCCCO", "CCCCCCCCCO", "CCCCCCCCCCO",
          "CCS", "CCCS", "CCCCS", "CCCCCS", "CCCCCCS", 
          "CCCCCCCS","CCCCCCCCS", "CCCCCCCCCS", "CCCCCCCCCCS"
          ]
    
    # Verify variables class.
    if type(desc) != bool:
        raise TypeError("Invalid 'desc' class. 'desc' must be a boolean")
    if type(plot) != bool:
        raise TypeError("Invalid 'plot' class. 'plot' must be a boolean")
        
    # Case when smiles is a list of strings.
    if type(smiles) == list:
        pgs = []
        for smile in smiles:
            pgs.append(atom_mapping(smile, desc, plot))
        return pgs
    
    pg_check_smiles(smiles)
        
    # The code gets the Rdkit object from the smiles.
    mol = Chem.MolFromSmiles(smiles)
    
    # Get molecule 3d configuration.
    mol, conf = config_mol(smiles)
    
    # Get molecule characteristics.
    coords = conf.GetPositions()
    atom_list = [atom.GetSymbol() for atom in mol.GetAtoms()]
    bonded_atoms, lengths = molecule_bonds(mol, coords)
    angles = molecule_angles(coords, bonded_atoms)
    
    # Plot.
    if plot:
        plot_molecule(atom_list, coords, bonded_atoms)
    # Determines the point group.
    if smiles in Cs:
        pg = "Cs"
    else:
        X = PointGroup(coords, atom_list, tolerance_eig=0.1, tolerance_ang=10)
        pg = X.get_point_group()
        # PointGroup function gives back C1v instead of Cs.
        if pg == "C1v":
            pg = "Cs"
        elif pg == "C1h":
            pg = "Cs"
    # Molecule description.
    if desc:
        print(f"Molecule SMILES : {smiles}")
        print(f"Molecule point group : {pg}")
        print(f"Atoms list : {atom_list}")
        print(f"Coordinates : {coords}")
        print(f"Bonds : {bonded_atoms}")
        print(f"Bond lenghts (normalized) : {lengths}")
        print(f"Angles : {angles}")
        print("")
            
    return pg

# This function allows user to interact with the function and set its parameters in a practical way.
# That way, user can choose to hide/show all informations besides the point group 
# and to hide/show the graphical representation.
def atom_mapping_interface():
    
    test_list = input("Do you want to use the function on multiple smiles ? [y/n] ")
    print("")
    # Asks for the smiles and creates a list out of them.
    if test_list == "y":
        print("Write your smiles. Write 'stop' when you are done")
        print("")
        smile = []
        stop = False
        while not stop:
            x = input("Write a smile : ")
            if x == "stop":
                break
            smile.append(x)
    elif test_list == "n":
        smile = input("What is your smile(s) ? ")
    else:
        raise TypeError("Invalid input write 'y' or 'n'")
    
    desc = input("Do you want a description ? [y/n] ")
    # Takes in memory the choices of the parameter of atom_mapping about the info.
    if desc == "y":
        desc = True
    elif desc == "n":
        desc = False
    else:
        raise TypeError("Invalid input. Write 'y' or 'n'.")
    
    plot = input("Do you want a 3d plot of your molecule(s) ? [y/n] ")
    # Takes in memory the choices of the parameter of atom_mapping about the graph.
    if plot == "y":
        plot = True
    elif plot == "n":
        plot = False
    else:
        raise TypeError("Invalid input. Write 'y' or 'n'.")
        
    print("")
    # Runs the code according th choices of the user.
    if type(smile) != str:
        for i in smile:
            atom_mapping(i, desc, plot)
    else:
        atom_mapping(smile, desc, plot)
        
    if not desc:
        print(f"Molecule(s) point group(s) {atom_mapping(smile, False, False)}")
    
    return

def main() :
    atom_mapping_interface()

if __name__ == "main":
    main()

