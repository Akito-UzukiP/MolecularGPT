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
        return pic
    except:
        print(smiles)
        raise Exception("Please ensure that the SMILES string is valid.")


def get_compound_properties(formula):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{formula}/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        return None


if __name__ == "__main__":
    formula = "C6H6"
    data = get_compound_properties(formula)

    if data is not None:
        print(json.dumps(data, indent=4))
    else:
        print("No data found.")
