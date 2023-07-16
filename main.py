import openai
import rdkit.Chem.Draw
import rdkit.Chem as Chem
import json
import os
import dash_bio as dashbio
import requests
from utils import _get_compound_properties, smiles_to_3d, cal_mol_props, pregenerated_molecule, docking_score, get_knowledge
from dash import Dash, dcc, html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import plotly.express as px
from dash.exceptions import PreventUpdate
import json
from dash import callback_context as ctx
from dash_bio.utils import xyz_reader
import dash_bio.utils.ngl_parser as ngl_parser
openai_api_key = os.environ.get("OPENAI_API_KEY")
openai.api_key = openai_api_key
system_prompt = open("system_prompt.txt","r",encoding='UTF-8').read()
main_history = []
welcome_word = open("./assets/welcome_word.txt","r",encoding='UTF-8').read()

init_smiles = 'CC(C1=CC2=C(C=C1)C3=CC(=C(C=C3OC2=O)OCC4=CC=C(C=C4)OC)OC)O'

# Keep the functions as they are...
functions = [
    {   
        "name": "pregenerated_molecule",
        "description": "使用深度学习模型生成一个高QED、符合靶点要求的分子,一次性生成1000个分子,从中选出药物性质评分最高的分子,如果用户要求生成分子则调用这个函数。请注意，这个分子会自动被其它组件调用生成图片，你不需要生成图片。",
        "parameters": {
            "type": "object",
            "properties": {
            },
        "required": []
        }
    },
    {
        "name": "_get_compound_properties",
        "description": "通过PubChem查询分子的属性。如果使用分子式查询，你必须用英语分子式子，,不是分子名字！用中文查询是非法的！如果通过SMILES查询只有分子式、SMILES和InChI表示，则你需要跟用户说明这个分子并没有收录，可能是全新的分子。请注意：使用SMILES式进行的查询比使用formula进行的精准的多得多得多！！如果用户已经提供了SMILES，请一定使用SMILES进行查询！",
        "parameters": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "The molecular's SMILES or formula."
                },
                "query_type": {
                    "type": "string",
                    "enum": ['formula', "smiles"],
                    "description": "Choose to search with SMILES or formula."
                }
            },
            "required": ["query", "query_type"]
        }
    },
        {
        "name": "cal_mol_props",
        "description": "Calculates and returns multiple properties of the molecule given its SMILES string. Uses the RDKit library.",
        "parameters": {
            "type": "object",
            "properties": {
                "smi": {
                    "type": "string",
                    "description": "The SMILES string of the molecule to calculate properties for. Please ensure that the SMILES string is valid."
                },
                "verbose": {
                    "type": "boolean",
                    "description": "A flag that, when set to true, makes the function print the calculated properties. Default value is false."
                }
            },
            "required": ["smi"]
        }
    },
    {   
        "name": "docking_score",
        "description": "使用深度学习计算药物分子和蛋白质靶点的结合评分。这个分数是通过AutoDock软件将分子和靶点进行对接计算得到的。",
        "parameters": {
            "type": "object",
            "properties": {
            },
        "required": []
        }
    },
    {
        "name": "get_knowledge",
        "description": "查询一些制药相关的知识，以及和本公司相关的信息。",
        "parameters": {
            "type": "object",
            "properties": {
                "key": {
                    "type": "string",
                    "description": "The key of the knowledge to query.",
                    "enum": ['靶点', '靶点配对', '分子设计', '分子合成', '临床试验', '靶点口袋']
                }
            },
            "required": ["key"]
        }
    }
    ]

global_info = {
    "history": [],  
    "system_prompt": open("system_prompt.txt","r",encoding='UTF-8').read(),
    "functions": functions,
    "input_text": "",
    "output_text": "",
    "mol_image": "",
    "displaying_molecule": init_smiles,
    "molecule_info": {},
    "molecule_evaluate": "分子量（MW）：406.4 氢键供体数（HBD）：1 氢键受体数（HBA）：6 LogP值：4.6 可旋转键数（RotB）：6 QED分数：0.37 手性中心数：1 极性表面积（TPSA）：78.1 BertzCT：1254.91",
    "function_response": "",

}


def chatbot(global_info):
    #从global_info中获取信息
    input_text = global_info["input_text"]
    history = global_info["history"]
    functions = global_info["functions"]
    system_prompt = global_info["system_prompt"]
    SMILES = global_info["displaying_molecule"]
    function_response = None
    #组建对话历史
    history.append({"role": "user", "content": input_text})
    temp_history = [{"role": "system", "content": system_prompt + "现在显示的分子的信息:" + str(global_info["molecule_info"]) + "现在显示的分子的SMILES:" + str(global_info['displaying_molecule'])}]
    for i in history:
        temp_history.append(i)
    #获取第一轮回复
    completion = openai.ChatCompletion.create(
        model="gpt-3.5-turbo-16k",
        messages=temp_history,
        functions=functions,
        temperature = 1)
    response_message = completion['choices'][0]['message']
    #print(response_message)
    #如果有函数调用，调用函数，进行第二轮对话
    if response_message.get("function_call"):
        available_functions = {
            "pregenerated_molecule": pregenerated_molecule,
            "_get_compound_properties": _get_compound_properties,
            "cal_mol_props": cal_mol_props,
            "docking_score": docking_score,
            "get_knowledge": get_knowledge,
        }
        function_name = response_message["function_call"]["name"]
        fuction_to_call = available_functions[function_name]
        function_args = json.loads(response_message["function_call"]["arguments"])
        print(function_name)
        print(function_args)
        function_response = fuction_to_call(
            **function_args
        )
        #如果function是get_compound_properties，把分子的信息存到global_info
        if function_name == "_get_compound_properties":
            global_info["molecule_info"] = function_response
        #如果function是生成分子，把分子的信息存到global_info
        if function_name == "pregenerated_molecule":
            global_info["displaying_molecule"] = function_response
        if function_name == "cal_mol_props":
            global_info["molecule_evaluate"] = function_response[-1]
            function_response = function_response[-1]
        temp_history.append(response_message)  
        temp_history.append(
            {
                "role": "function",
                "name": function_name,
                "content": str(function_response),
            }
        )
        second_response = openai.ChatCompletion.create(
            model="gpt-3.5-turbo-16k",
            messages=temp_history,
            temperature = 1)
        second_response_message = second_response['choices'][0]['message']['content']
        history.append({"role": "assistant", "content": second_response_message})
    else:
        history.append({"role": "assistant", "content": response_message["content"]})

    #解析输出，根据AI的指令调用函数

    #将对话的结果返回给global_info
    global_info["history"] = history
    global_info["function_response"] = function_response
    global_info["output_text"] = response_message["content"]
    #print(global_info)
    return global_info

# Create a Dash app



app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP],title='分子可视化平台')
#CORS(app.server)  # 允许所有源的请求
# Define the layout
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H2('MolecularGPT Chatbot', className='text-center text-primary mb-4'),
        ], width=12, className='mx-auto')
    ], className='mx-auto'),
    dbc.Row([
        #第一列（chatGPT交流、分子编辑
        dbc.Col([
            html.P('当你在JSME上画好分子之后，你就可以右键->Copy as SMILES -> 复制分子信息 -> 粘贴到画面中间的分子结构表示法(SMILES)上，然后你就可以看到分子的模型啦！', style={'color': 'green' , 'font-size': '10px'}),
            dashbio.Jsme(id = 'jsme',width='100%',height='30vh',options="NOexportInChI NOexportInChIauxInfo NOexportInChIkey NOsearchInChIkey NOexportSVG paste newLook NOshowDragAndDrop SymbolInDepictMode", smiles=init_smiles, style={'width': '50vh', 'height': '30vh', 'border': 'solid grey 3px'}),
            dbc.Card(id='chat-history', children=[
                dbc.Card([
                    dbc.CardBody([
                        html.Img(src='./assets/chatgpt.png', className='icon-class'),  # You should replace /path/to/icon.png with the actual path to your icon image
                        dcc.Markdown(welcome_word, className='card-text',style={'color': 'white'}),
                    ], className='card bg-primary text-white mb-2', style={'display': 'flex','flex-direction':'row', 'width':'100%',  'align-items': 'center'})
                ], className='mb-4')
            ], className='mt-4', style={ 'align-items': 'center', 'height': '40vh', 'width':'50vh', 'border': 'solid grey 3px', 'overflowY': 'scroll','back-ground-color':'grey'}),
            dbc.FormGroup([
                dbc.Label('用户输入', className='form-label' , style={'color': 'white'}),
                dbc.Textarea(id='input-text', value="你可以为我做什么？", className='form-control', style={'align-items': 'center', 'border': 'solid grey 3px','width':'50vh', 'height': '5vh'})
                ]),
            dbc.Row([
            dbc.Button('发送内容', id='submit-button', n_clicks=0, color='primary', className='mt-2'),
            dbc.Button('生成分子', id='generation-button', n_clicks=0, color='primary', className='mt-2'),
            dbc.Button('查询分子', id='query-button', n_clicks=0, color='primary', className='mt-2'),
            dbc.Button('评价分子', id='eval-button', n_clicks=0, color='primary', className='mt-2'),
            dbc.Button('靶点结合评分', id='dock-eval-button', n_clicks=0, color='primary', className='mt-2'),
        ]),
        ], width=4),
        #第二列（Speck显示、分子属性显示
        dbc.Col([
            dbc.Row([
                dashbio.Speck(
                    id='mol-Speck',
                    data=smiles_to_3d(init_smiles),
                    view={
                        'resolution': 500,
                        'ao': 0.1,
                        'outline': 1,
                        'atomScale': 0.25,
                        'relativeAtomScale': 0.33,
                        'bonds': True,
                        'zoom' : 0.05
                    },
                    style={
                        'height': '500px',
                        'width': '500px'
                    },
                )],
                style={
                    'background-image': './assets/PRDM9.png',
                    'background-size': 'cover',
                }),
            dbc.Row([
                dbc.Col([
                    dbc.Label("分子表示方式", className="form-label", style={'margin-top': '10px', 'color': 'white'}),
                    dcc.Dropdown(
                        id='speck-style-dropdown',
                        options=[
                            {'label': '球模型', 'value': 'default'},
                            {'label': '球棍模型', 'value': 'stickball'},
                            {'label': '卡通效果', 'value': 'toon'},
                            {'label': '棍棒模型', 'value': 'licorice'},
                            # Add more styles here
                        ],
                        value='toon',  # Default value
                        style={'width': '50vh', 'margin-top': '10px'}
                    ),
                    dbc.FormGroup([
                        dbc.Label("分子结构表示法(SMILES)", className="form-label", style={'margin-top': '10px', 'color': 'white'}),
                        dbc.Textarea(
                            id="smiles-textbox",
                            value=init_smiles,
                            style={'height': '5vh', 'color': 'blue'},
                            readOnly=False
                        )
                    ], style={'width': '50vh', 'margin-top': '10px'}),

                    dbc.FormGroup([
                        dbc.Label("分子评分", className="form-label", style={'margin-top': '10px', 'color': 'white'}),
                        dcc.Markdown(
                            id="molecule-properties",
                            children="分子量（MW）：406.4 氢键供体数（HBD）：1 氢键受体数（HBA）：6 LogP值：4.6 可旋转键数（RotB）：6 QED分数：0.37 手性中心数：1 极性表面积（TPSA）：78.1 BertzCT：1254.91",
                            style={'color': 'white'}
                        )
                    ],style={'width': '50vh', 'margin-top': '10px'})
                ], width=4)
            ], className='mx-auto'),

        ], width=4),
        #第三列
        dbc.Col([
            html.H2('靶点蛋白质SH3b', className='text-center text-primary mb-4',style={'border': 'solid grey 3px', 'background-color': 'white'}),
            html.Img(id = 'target_png', src='./assets/target.png',style={'width': '500px', 'height': '500px'}),
            dbc.Button('显示真实药物分子', id='show-answer', n_clicks=0, color='primary', className='mt-2')
        ], width=4, className='mx-auto'),
    ], className='mx-auto'),

    dcc.Store(id='global-store',data=global_info),
    dcc.Store(id='temp-store',data=global_info)
], fluid=True, className='mx-auto')






# chatGPT输入（文本框输入，聊天栏输出，global_store存储）
@app.callback(
    Output('chat-history', 'children'),
    Output('smiles-textbox', 'value'),
    Output('temp-store', 'data'),
    Input('submit-button', 'n_clicks'),
    Input('query-button', 'n_clicks'),
    Input('generation-button', 'n_clicks'),
    Input('eval-button', 'n_clicks'),
    Input('dock-eval-button', 'n_clicks'),
    State('input-text', 'value'),
    State('global-store', 'data')
)
def update_chat(n_clicks,n_click_query,n_click_generation,n_click_eval,n_click_dock_eval, input_text, global_store):
    if not ctx.triggered:
        raise PreventUpdate
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # 根据不同的按钮设置不同的 input_text
    print(button_id)
    if button_id == 'query-button':
        input_text = "请帮我用PubChem查询一下现在这个分子的信息"
    elif button_id == 'generation-button':
        input_text = "请用深度学习模型生成一个高评分、符合靶点要求的分子"
    elif button_id == 'eval-button':
        input_text = "请评价一下这个分子"
    elif button_id == 'dock-eval-button':
        input_text = "请评价一下目前这个分子的靶点结合评分如何"

    if n_clicks > 0 or n_click_query > 0 or n_click_generation > 0 or n_click_eval > 0 or n_click_dock_eval > 0:
        
        global_store["input_text"] = input_text
        #更新global_store
        datas = chatbot(global_store)
        chat_history = datas['history']
        chat_output = []
        user_input = []
        bot_output = []
        for message in chat_history:
            if message['role'] == 'user':
                user_input.append(message['content'])
            else:
                bot_output.append(message['content'])
        #print(user_input, bot_output)
        for i in range(len(user_input)):
            chat_output.append(
                dbc.Card([
                    dbc.CardBody([
                        html.Img(src='./assets/jingtai.png', className='icon-class'),  # You should replace /path/to/icon.png with the actual path to your icon image
                        dcc.Markdown(user_input[i], className='card-text',style={'color': 'grey'}),
                    ], className='card bg-light text-dark mb-2', style={'display': 'flex', 'flex-direction':'row', 'width':'95%',  'align-items': 'center', 'word-wrap': 'break-word', 'overflow-wrap': 'break-word'}),
                    dbc.CardBody([
                        html.Img(src='./assets/chatgpt.png', className='icon-class'),  # You should replace /path/to/icon.png with the actual path to your icon image
                        dcc.Markdown(bot_output[i], className='card-text',style={'color': 'white'}),
                    ], className='card bg-primary text-white mb-2', style={'display': 'flex','flex-direction':'row', 'width':'95%',  'align-items': 'center', 'word-wrap': 'break-word', 'overflow-wrap': 'break-word'})
                ], className='mb-4')
            )

        #print(chat_output)
        return chat_output, datas['displaying_molecule'], datas
    return [],"",  global_store


#textbox传输到mol-Speck

@app.callback(
    Output('mol-Speck', 'data'),
    Input('global-store', 'data'), 
    State('mol-Speck', 'data')
)
def update_mol_speck(global_store,current_data):
    smiles = global_store["displaying_molecule"]
    if smiles is not None and smiles != "":
        #确认SMILES是否合法
        temp = smiles_to_3d(smiles)
        if temp is not None:
            print(smiles)
            return temp
    return current_data

#textbox传输到ngl
'''@app.callback(
    Output('default-ngl-molecule', 'data'),
    Input('global-store', 'data'), 
    State('default-ngl-molecule', 'data')
)
def update_mol_ngl(global_store,current_data):
    smiles = global_store["displaying_molecule"]
    if smiles is not None and smiles != "":
        #确认SMILES是否合法
        temp = smiles2pdb(smiles)
        current_data[0]['config']['input'] = temp

    return current_data'''


#全局信息-->jsme/textbox显示

@app.callback(
    Output('jsme', 'smiles'),
    Input('global-store', 'data')  # The data of mol-Speck is updated when the global store is updated
)
def update_jsme(global_store):
    return Chem.MolToSmiles(Chem.MolFromSmiles(global_store["displaying_molecule"]))

#style变更

@app.callback(
    Output('mol-Speck', 'presetView'),
    Input('speck-style-dropdown', 'value'),
)
def update_speck_style(style_value):
    # Define different styles
    styles = {
        'default': 'default',
        'stickball': 'stickball',
        'toon': 'toon',
        'licorice': 'licorice',
        # Add more styles here
    }

    return styles[style_value]

#显示分子属性

@app.callback(
    Output('molecule-properties', 'children'),
    Input('global-store', 'data')
)
def update_molecule_properties(global_store):
    return global_store["molecule_evaluate"]

#从smiles-textbox更新global-store
@app.callback(
    Output('global-store', 'data'),
    Input('smiles-textbox', 'value'),
    State('temp-store', 'data')
)
def update_global_store_from_smiles_textbox(smiles, global_store):
    global_store["displaying_molecule"] = smiles
    global_store['molecule_evaluate'] = cal_mol_props(smiles)[-1]
    return global_store

@app.callback(
    Output('target_png', 'src'),
    Input('show-answer', 'n_clicks')
)
def update_target_png(n_clicks):
    if n_clicks % 2 == 0:
        return './assets/target.png'
    else:
        return './assets/target+drug.png'
# Run the app
if __name__ == '__main__':
    #serve(app.server,port=8050)
    app.run_server(debug=True, host='0.0.0.0', port=8050)
