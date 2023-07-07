import openai
import rdkit.Chem.Draw
import rdkit.Chem as Chem
import gradio as gr
import json
import os
import requests
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
# ["Compound", "Compound Complexity", "Count", "Fingerprint", "IUPAC Name", "InChI", "InChiKey", "Log P", "Mass", "Molecular Formula", "Molecular Weight", "SMILES", "TOPOLOGICAL", "Weight", "ALL"]
def get_compound_properties(query, query_type='formula',prop_type='Molecular Weight'):
    props = _get_compound_properties(query, query_type)
    if props is None:
        return "No data found."
    if prop_type == "ALL":
        return props
    return_props = []
    for mole in props:
        prop_dicts = {}
        for i in mole:
            if i == prop_type:
                prop_dicts[i] = mole[i]
        return_props.append(prop_dicts)
    return return_props
    




if __name__ == "__main__":
    formula = "C6H6"
    data = get_compound_properties(formula)

