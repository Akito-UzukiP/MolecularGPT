# MolecularGPT
## 使用GPT-3.5 结合各个工具的分子交互软件
### 1. 介绍
![Alt text](image.png)
### 2. 安装
- 克隆本repo后直接使用pip安装即可, openbabel需要使用conda安装。注意有些package需要MSVC编译器，如果没有的话会报错。
```bash
git clone https://github.com/Akito-UzukiP/MolecularGPT.git
cd MolecularGPT
pip install -r requirements.txt

conda install -c conda-forge openbabel 
```
- 请注意，必须正确设置系统conda环境变量，否则openbabel无法正常运行。

- 需要设置openai的api key环境变量，或者直接在代码里设置，以下代码是通过PowerShell命令行设置API_KEY，将sk-******更换成你的API KEY
```PowerShell
[Environment]::SetEnvironmentVariable("OPENAI_API_KEY", "sk-******", "User")
```

### 3. 使用
目前实现的自然语言调用功能：
- 调用GCPN预训练模型/GCPN疏水强化学习模型生成分子（为方便展示，预先生成了1000个分子以供选取，相关代码在generator.py中，模型权重未包括[其实是不小心删了])
- 通过SMILES或者分子式询问PubChem数据库，获取分子的物化性质（同分子式的分子可能有多个异构体，只选取第一个，建议使用SMILES查询）
- 给予重要关键词信息库，在chatGPT中询问时它会自动去调用。
- 使用Dash-Bio完成画面设计。

## TODO:
- 实现真实的AutoDock4分子对接计算。
- 实现真实的靶点显示。
