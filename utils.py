import openai
import rdkit.Chem.Draw
import rdkit.Chem.AllChem as Chem
from openbabel import openbabel, pybel
import gradio as gr
import json
import os
import requests
from dash_bio.utils import xyz_reader
from generator import GCPN_simple_molecule_generation,GCPN_hydrophobic_molecule_generation

def display_mol(smiles):
    #全大写
    smiles = smiles.upper()
    try:
        mol = Chem.MolFromSmiles(smiles)
        pic = rdkit.Chem.Draw.MolToImage(mol)
        print(pic)
        pic.save("2D_pic.png")
        return pic
    except:
        print(smiles)
        raise Exception("Please ensure that the SMILES string is valid.")



def smiles_to_formula(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    return Chem.rdMolDescriptors.CalcMolFormula(molecule)

# 获取pubchem的分子属性
def __get_compound_properties(query, query_type='formula'):
    if query_type.lower() == 'formula':
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{query}/JSON"
    elif query_type.lower() == 'smiles':
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{query}/JSON"
    else:
        print(f"Invalid query type: {query_type}. Please choose 'formula' or 'smiles'.")
        return None

    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        return data
    else:
        return None

def _get_compound_properties(query, query_type='formula'):
    data = __get_compound_properties(query, query_type)
    if data is  None:
        return None
    props = []
    for mole in data['PC_Compounds']:
        prop_dict = {}
        if mole['id'].get('id') is not None:
            prop_dict['cid'] = mole['id']['id']['cid']
        for prop in data['PC_Compounds'][0]['props']:
            prop_dict[prop['urn']['label']] = prop['value']
        props.append(prop_dict)
    return props
    
#支持的prop_type:
# Compound
# Compound Complexity
# Count
# Fingerprint
# IUPAC Name
# InChI
# InChIKey
# Log P
# Mass
# Molecular Formula
# Molecular Weight
# SMILES
# Topological
# Weight
# ["Compound", "Compound Complexity", "Count", "Fingerprint", "IUPAC Name", "InChI", "InChiKey", "Log P", "Mass", "Molecular Formula", "Molecular Weight", "SMILES", "TOPOLOGICAL", "Weight"]
def get_compound_properties(query, query_type='formula',prop_type=['SMILES']):
    props = _get_compound_properties(query, query_type)
    if props is None:
        return "No data found."
    return_props = []
    for mole in props:
        prop_dicts = {}
        for i in mole:
            for j in prop_type:
                if i == j:
                    prop_dicts[i] = mole[i]
        return_props.append(prop_dicts)
    return return_props

def interact_with_chatgpt(messages, functions, model="gpt-3.5-turbo-0613",temperature=0.5):
    openai.api_key = os.getenv("OPENAI_API_KEY")
    response = openai.ChatCompletion.create(
        model=model,
        messages=messages,
        functions=functions,
        temperature = temperature)
    return response['choices'][0]['message']


def smiles_to_3d(smiles):
    # Convert the SMILES string to a RDKit molecule
    mol = Chem.MolFromSmiles(smiles)

    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Generate a 3D conformer for the molecule
    Chem.EmbedMolecule(mol)
    
    # Convert the RDKit molecule to a Open Babel molecule
    obmol = openbabel.OBMol()
    converter = openbabel.OBConversion()
    converter.SetInAndOutFormats('mol', 'xyz')
    converter.ReadString(obmol, Chem.MolToMolBlock(mol))

    # Write the Open Babel molecule to a .xyz file
    #with open(filename, 'w') as f:
    #    f.write(pybel.Molecule(obmol).write('xyz'))
    return xyz_reader.read_xyz(pybel.Molecule(obmol).write('xyz'), is_datafile=False)

def combine_molecules(smiles1, smiles2, atom_index1, atom_index2):
    # 从 SMILES 字符串创建分子
    molecule1 = Chem.MolFromSmiles(smiles1)
    molecule2 = Chem.MolFromSmiles(smiles2)

    # 创建一个反应，用于将两个分子在指定的原子上合并
    reaction = Chem.ReactionFromSmarts('[*:1].[*:2]>>([*:1].[*:2])')

    # 执行反应
    product = reaction.RunReactants((molecule1, molecule2))

    # 从产物中获取新的分子
    new_molecule = product[0][0]

    # 返回新分子的 SMILES 字符串
    return Chem.MolToSmiles(new_molecule)


if __name__ == "__main__":
    # 甲烷和乙烷的 SMILES 字符串
    smiles1 = 'C'
    smiles2 = 'CC'

    # 合并甲烷和乙烷，得到丙烷
    new_smiles = combine_molecules(smiles1, smiles2, 0, 0)

    print(new_smiles)  # 输出：'CCC'

