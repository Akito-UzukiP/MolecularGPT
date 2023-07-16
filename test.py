from dash_bio.utils import xyz_reader
from openbabel import pybel
import openbabel
import rdkit.Chem.AllChem as Chem
from rdkit.Chem import Draw
def smiles_to_3d(smiles):
    # Convert the SMILES string to a RDKit molecule
        print(smiles)
        mol = Chem.MolFromSmiles(smiles)
        print(mol)
        # Add hydrogens to the molecule
        mol = Chem.AddHs(mol)
        print(mol)
        # Generate a 3D conformer for the molecule
        Chem.EmbedMolecule(mol)
        print(mol)       
        # 显示mol
        Draw.MolToImage(mol,show=True)
        # Convert the RDKit molecule to a Open Babel molecule
        obmol = openbabel.OBMol()
        converter = openbabel.OBConversion()
        print(converter)
        converter.SetInAndOutFormats('mol', 'xyz')
        converter.ReadString(obmol, Chem.MolToMolBlock(mol))
        print(pybel.Molecule(obmol).write('xyz'))
        # Write the Open Babel molecule to a .xyz file
        #with open(filename, 'w') as f:
        #    f.write(pybel.Molecule(obmol).write('xyz'))
        return xyz_reader.read_xyz(pybel.Molecule(obmol).write('xyz'), is_datafile=False)
smiles_to_3d('c1ccccc1')