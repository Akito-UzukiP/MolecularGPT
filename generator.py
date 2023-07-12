import torch
from torchdrug import datasets, core, models, tasks, utils
from torch import optim
import warnings
warnings.filterwarnings("ignore")
import logging
logging.disable(logging.CRITICAL)
import pickle
from utils import cal_mol_props
import numpy as np
from rdkit import Chem
# 用GCPN生成随机符合规则的分子
def GCPN_simple_molecule_generation(num=1):
    with open("./torchdrug/zinc250k.pkl", "rb") as fin:
        dataset = pickle.load(fin)
    dataset.transform = None
    model = models.RGCN(input_dim=18,
                        num_relation=3,
                        hidden_dims=[256, 256, 256, 256], batch_norm=False)
    task = tasks.GCPNGeneration(model, [6, 7, 8, 9, 15, 16, 17, 35, 53], max_edge_unroll=12,
                                max_node=38, criterion="nll")
    optimizer = optim.Adam(task.parameters(), lr = 1e-3)
    solver = core.Engine(task, dataset, None, None, optimizer,
                             batch_size=128, log_interval=10)
    solver.load("./torchdrug/gcpn_zinc250k_qed_rl_1epoch.pkl")
    results = task.generate(num_sample=num, max_resample=10)
    return results.to_smiles()

# 用GCPN生成疏水性良好的的分子
def GCPN_hydrophobic_molecule_generation(num=1):
    with open("./torchdrug/zinc250k.pkl", "rb") as fin:
        dataset = pickle.load(fin)
    dataset.transform = None
    model = models.RGCN(input_dim=18,
                        num_relation=3,
                        hidden_dims=[256, 256, 256, 256], batch_norm=False)
    task = tasks.GCPNGeneration(model, [6, 7, 8, 9, 15, 16, 17, 35, 53], max_edge_unroll=12,
                                max_node=38, criterion="nll")
    optimizer = optim.Adam(task.parameters(), lr = 1e-3)
    solver = core.Engine(task, dataset, None, None, optimizer,
                             batch_size=128, log_interval=10)
    solver.load("./torchdrug/gcpn_zinc250k_rl_1epoch.pkl")
    results = task.generate(num_sample=num, max_resample=10)
    return results.to_smiles()

def generate_good_molecule():
    # 用GCPN_simple_molecule_generation生成20个，取其中QED最好的。
    moles = GCPN_simple_molecule_generation(100)
    mole_qed = []
    for i in moles:
        mole_qed.append(cal_mol_props(i, verbose=False)[3])
    #排序，找到top-1的位置
    mole_qed = np.array(mole_qed)
    mole_qed = np.argsort(mole_qed)
    return Chem.MolToSmiles(Chem.MolFromSmiles(moles[mole_qed[-1]]))

