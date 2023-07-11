import openai
import rdkit.Chem.Draw
import rdkit.Chem.AllChem as Chem
from openbabel import openbabel, pybel
import json
import os
import requests
from dash_bio.utils import xyz_reader
import numpy as np
from rdkit.Chem import QED, Descriptors, rdMolDescriptors
from scipy.stats import gaussian_kde
import rdkit
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
    return props[0]
    
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
'''def get_compound_properties(query, query_type='formula',prop_type=['SMILES']):
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
    return return_props'''

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
    try:
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
    except:
        return None
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

def cal_mol_props(smi, verbose=False):
    try:
        m = Chem.MolFromSmiles(smi)
        if not m:
            return None, None, None, None, None, None, None, None, None
 
        logp = np.round(Descriptors.MolLogP(m), 2)
        tpsa = np.round(Descriptors.TPSA(m), 1)
        mw = np.round(Descriptors.MolWt(m), 1)
        qed = np.round(QED.qed(m), 2)
        hba = rdMolDescriptors.CalcNumLipinskiHBA(m)
        hbd = rdMolDescriptors.CalcNumLipinskiHBD(m)
        rob = rdMolDescriptors.CalcNumRotatableBonds(m)
        chiral_center = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
 
        # 计算Bertz CT的数据分布的直方图
        bertz_ct = Descriptors.BertzCT(m)
 
 
        if verbose:
            print(smi)
            print('MW ', mw)
            print('HBD ', hbd)
            print('HBA ', hba)
            print('Logp ', logp)
            print('RotB ', rob)
            print('QED ', qed)
            print('chiral_center ', chiral_center)
            print('TPSA ', tpsa)
            print('bertz_ct', bertz_ct)
 
        return logp, tpsa, mw, qed, hba, hbd, rob, chiral_center, bertz_ct,"分子的SMILES表示为：{}，其各项性质如下：\n\
分子量（MW）：{}\n\
氢键供体数（HBD）：{}\n\
氢键受体数（HBA）：{}\n\
LogP值：{}\n\
可旋转键数（RotB）：{}\n\
QED分数：{}\n\
手性中心数：{}\n\
极性表面积（TPSA）：{}\n\
BertzCT：{}".format(smi, mw, hbd, hba, logp, rob, qed, chiral_center, tpsa, bertz_ct)
 
    except Exception as e:
        print(e)
        return None, None, None, None, None, None, None, None, None, None


if __name__ == "__main__":
    # 甲烷和乙烷的 SMILES 字符串
    SMILES = 'CCCc1ccc(Cc2sc3c(c2C(=O)NC(C)c2ccc(C(=O)O)cc2)CCOC3)cc1'
    logp, tpsa, mw, qed, hba, hbd, rob, chiral_center, bertz_ct, word = cal_mol_props(SMILES, verbose=False)
    print(word)

