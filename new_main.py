import openai
import rdkit.Chem.Draw
import rdkit.Chem as Chem
import json
import os
import dash_bio as dashbio
import requests
from generator import GCPN_hydrophobic_molecule_generation, generate_good_molecule
from utils import _get_compound_properties, display_mol, smiles_to_3d, cal_mol_props
from dash import Dash, dcc, html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import plotly.express as px
import json
from dash_bio.utils import xyz_reader
openai_api_key = os.environ.get("OPENAI_API_KEY")
openai.api_key = openai_api_key
system_prompt = open("system_prompt.txt","r",encoding='UTF-8').read()
main_history = []
proxies = {
        "http":  "socks5h://localhost:11796",
        "https": "socks5h://localhost:11796",
}
requests.Session().proxies = proxies




# Keep the functions as they are...
functions = [
    {
        "name": "display_mol",
        "description": "显示分子结构，输入分子正确的SMILES字符串，通过RDKit库绘制分子结构2D图像,直接显示在外面，不会返回给AI",
        "parameters": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "string",
                    "description": "The SMILES string of the molecule to display. Please ensure that the SMILES string is valid."
                }
            },
        "required": ["smiles"]
        }
    },
    {   
        "name": "generate_good_molecule",
        "description": "生成一个高QED的分子",
        "parameters": {
            "type": "object",
            "properties": {
            },
        "required": []
        }
    },
    {
        "name": "_get_compound_properties",
        "description": "Get compound information from PubChem. Choose which property to get. If you get None, it means that this compond is a de novo generated compound.",
        "parameters": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "The molecular's SMILES or formula."
                },
                "query_type": {
                    "type": "string",
                    "enum": ["formula", "smiles"],
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
    }
    ]

global_info = {
    "history": [],  
    "system_prompt": open("system_prompt.txt","r",encoding='UTF-8').read(),
    "functions": functions,
    "input_text": "",
    "output_text": "",
    "mol_image": "",
    "xyz_file": {},
    "molecule_info": {},
    "molecule_evaluate": "",
    "function_response": ""
}


def chatbot(global_info):
    #从global_info中获取信息
    input_text = global_info["input_text"]
    history = global_info["history"]
    functions = global_info["functions"]
    system_prompt = global_info["system_prompt"]
    function_response = None
    #组建对话历史
    history.append({"role": "user", "content": input_text})
    temp_history = [{"role": "system", "content": system_prompt}]
    for i in history:
        temp_history.append(i)
    #获取第一轮回复
    completion = openai.ChatCompletion.create(
        model="gpt-3.5-turbo-0613",
        messages=temp_history,
        functions=functions,
        temperature = 1)
    response_message = completion['choices'][0]['message']
    #print(response_message)
    #如果有函数调用，调用函数，进行第二轮对话
    if response_message.get("function_call"):
        available_functions = {
            "generate_good_molecule": generate_good_molecule,
            "GCPN_hydrophobic_molecule_generation": GCPN_hydrophobic_molecule_generation,
            "_get_compound_properties": _get_compound_properties,
            "cal_mol_props": cal_mol_props
        }
        function_name = response_message["function_call"]["name"]
        fuction_to_call = available_functions[function_name]
        function_args = json.loads(response_message["function_call"]["arguments"])
        #print(function_args)
        function_response = fuction_to_call(
            **function_args
        )
        #如果function是get_compound_properties，把分子的信息存到global_info
        if function_name == "_get_compound_properties":
            global_info["molecule_info"] = function_response
        #如果function是生成分子，把分子的信息存到global_info
        if function_name == "generate_good_molecule" or function_name == "GCPN_hydrophobic_molecule_generation":
            global_info["xyz_file"] = {'SMILES':function_response}
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
            model="gpt-3.5-turbo-0613",
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




app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Define the layout
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H2('MolecularGPT Chatbot', className='text-center text-primary mb-4'),
        ], width=12)
    ]),
    dbc.Row([
        dbc.Col([
            dbc.Card(id='chat-history', children=[], className='mt-4', style={'align-items': 'center', 'height': '800px', 'width':'800px', 'border': 'solid grey 3px', 'overflowY': 'scroll'}),
            dbc.FormGroup([
                dbc.Label('用户输入', className='form-label'),
                dbc.Textarea(id='input-text', className='form-control', style={'align-items': 'center', 'border': 'solid grey 3px','width':'800px', 'height': '100px'}),
                dbc.FormText('请输入你的查询.', color='secondary'),
                ]),
            dbc.Button('提交', id='submit-button', n_clicks=0, color='primary', className='mt-2'),
        ], width=6),
        dbc.Col([
            dcc.Loading(
                id="loading",
                type="cube",
                children=dashbio.Speck(
                    id='mol-Speck',
                    data = xyz_reader.read_xyz('./assets/test.xyz'),
                    view={
                        'resolution': 800,
                        'ao': 0.1,
                        'outline': 1,
                        'atomScale': 0.25,
                        'relativeAtomScale': 0.33,
                        'bonds': True,
                        'zoom' : 0.05
                    },
                    style={
                        'height': '800px',
                        'width': '800px',
                        'border': 'solid grey 3px'
                    }
                )
            ),
            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='speck-style-dropdown',
                        options=[
                            {'label': 'Style 1', 'value': 'style1'},
                            {'label': 'Style 2', 'value': 'style2'},
                            {'label': 'Style 3', 'value': 'style3'},
                            # Add more styles here
                        ],
                        value='style1'  # Default value
                    ),
                    dbc.FormGroup([
                        dbc.Label("properties", className="form-label", style={'margin-top': '10px'}),
                        dbc.Textarea(
                            id="molecule-properties",
                            value="",
                            style={'height': '150px'},
                            readOnly=True,
                        )
                    ])
                ], width=6)
            ]),

        ], width=6)
    ]),
    dcc.Store(id='global-store',data=global_info)
], fluid=True)






# Define the callback
@app.callback(
    Output('chat-history', 'children'),
    Output('global-store', 'data'),
    Input('submit-button', 'n_clicks'),
    State('input-text', 'value'),
    State('global-store', 'data')
)
def update_chat(n_clicks, input_text, global_store):
    if n_clicks > 0:
        global_store["input_text"] = input_text
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
                        html.P("   "+user_input[i], className='card-text'),
                    ], className='card bg-light text-dark mb-2',style={'display': 'flex', 'flex-direction':'row',  'align-items': 'center'}),
                    dbc.CardBody([
                        html.Img(src='./assets/chatgpt.png', className='icon-class'),  # You should replace /path/to/icon.png with the actual path to your icon image
                        html.P("   "+bot_output[i], className='card-text'),
                    ], className='card bg-primary text-white mb-2',style={'display': 'flex','flex-direction':'row',  'align-items': 'center'})
                ], className='mb-4')
            )
        #print(chat_output)
        return chat_output, datas
    return [], global_store


@app.callback(
    Output('mol-Speck', 'data'),
    Input('global-store', 'data'),  # The data of mol-Speck is updated when the global store is updated
    State('submit-button', 'n_clicks'),
    State('mol-Speck', 'data')
)
def update_mol_speck(global_store, n_clicks,current_data):
    if n_clicks == 0:
        return current_data
    print(global_store["xyz_file"])
    if global_store["xyz_file"] is not None and global_store["xyz_file"].get("SMILES"):
        return smiles_to_3d(global_store["xyz_file"]["SMILES"])
    return current_data

@app.callback(
    Output('mol-Speck', 'view'),
    Input('speck-style-dropdown', 'value'),
)
def update_speck_style(style_value):
    # Define different styles
    styles = {
        'style1': {
            'resolution': 800,
            'ao': 0.1,
            'outline': 1,
            'atomScale': 0.25,
            'relativeAtomScale': 0.33,
            'bonds': True,
            'zoom' : 0.05
        },
        'style2': {
            'resolution': 800,
            'ao': 0.1,
            'outline': 0,
            'atomScale': 0.1,
            'relativeAtomScale': 0.7,
            'bonds': True,
            'zoom' : 0.05,
            'adomShade' : 1
        },
        'style3': {
            'resolution': 800,
            'ao': 0.1,
            'outline': 0,
            'atomScale': 0.8,
            'relativeAtomScale': 0.9,
            'bonds': False,
            'zoom' : 0.05
        }
        # Add more styles here
    }

    return styles[style_value]

@app.callback(
    Output('molecule-properties', 'value'),
    Input('global-store', 'data')
)
def update_molecule_properties(global_store):
    return global_store["molecule_evaluate"]

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)
