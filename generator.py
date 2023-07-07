import torch
from torchdrug import datasets, core, models, tasks, utils
from torch import optim
import warnings
warnings.filterwarnings("ignore")
import logging
logging.disable(logging.CRITICAL)
import pickle


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
                            gpus=(0,), batch_size=128, log_interval=10)
    solver.load("./torchdrug/gcpn_zinc250k_5epoch.pkl")
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
                            gpus=(0,), batch_size=128, log_interval=10)
    solver.load("./torchdrug/gcpn_zinc250k_rl_1epoch.pkl")
    results = task.generate(num_sample=num, max_resample=10)
    return results.to_smiles()

