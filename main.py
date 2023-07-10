import openai
import rdkit.Chem.Draw
import rdkit.Chem as Chem
import gradio as gr
import json
import os
import requests
from generator import GCPN_simple_molecule_generation,GCPN_hydrophobic_molecule_generation
from utils import get_compound_properties, display_mol
openai_api_key = os.environ.get("OPENAI_API_KEY")
openai.api_key = openai_api_key
system_prompt = open("system_prompt.txt","r",encoding='UTF-8').read()
main_history = []
proxies = {
        #          [协议]://  [地址]  :[端口]
        "http":  "socks5h://localhost:11796",  # 再例如  "http":  "http://127.0.0.1:7890",
        "https": "socks5h://localhost:11796",  # 再例如  "https": "http://127.0.0.1:7890",
}
requests.Session().proxies = proxies
#设定requests的代理



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
        "name": "GCPN_simple_molecule_generation",
        "description": "用GCPN生成随机符合规则的分子",
        "parameters": {
            "type": "object",
            "properties": {
                "num": {
                    "type": "integer",
                    "description": "生成分子的数量,默认为1"
                },
            },
        "required": []
        }
    },
    {
        "name": "GCPN_hydrophobic_molecule_generation",
        "description": "用GCPN生成随机符合规则的疏水分子",
        "parameters": {
            "type": "object",
            "properties": {
                "num": {
                    "type": "integer",
                    "description": "生成分子的数量,默认为1"
                },
            },
        },
        "required": []
    },
    {
        "name": "get_compound_properties",
        "description": "Get compound information from PubChem. Choose which property to get. If you want to get all roperties, use 'ALL'",
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
                },
                "prop_type":{
                    "type": "string",
                    "enum":  ["Compound", "Compound Complexity", "Count", "Fingerprint", "IUPAC Name", "InChI", "InChIKey", "Log P", "Mass", "Molecular Formula", "Molecular Weight", "SMILES", "Topological", "Weight", "ALL"],
                    "description": "Which property you want to get. Use 'ALL' to get all properties"
                }
            },
            "required": ["query", "query_type", "prop_type"]
        }
    }
    ]




def chatbot(input_text):
    function_response = None
    main_history.append({"role": "user", "content": input_text})
    temp_history = []
    temp_history = [{"role": "system", "content": system_prompt}]
    for i in main_history:
        temp_history.append(i)
    completion = openai.ChatCompletion.create(
        model="gpt-3.5-turbo-0613",
        messages=temp_history,
        functions=functions,
        temperature = 1)
    response_message = completion['choices'][0]['message']
    print(response_message)
    if response_message.get("function_call"):
        available_functions = {
            "display_mol": display_mol,
            "GCPN_simple_molecule_generation": GCPN_simple_molecule_generation,
            "GCPN_hydrophobic_molecule_generation": GCPN_hydrophobic_molecule_generation,
            "get_compound_properties": get_compound_properties
        }
        function_name = response_message["function_call"]["name"]
        fuction_to_call = available_functions[function_name]
        function_args = json.loads(response_message["function_call"]["arguments"])
        print(function_args)
        function_response = fuction_to_call(
            **function_args
        )
        if function_name == "display_mol":
            function_response.save("2D_pic.png")
            function_response = "PIC displayed"

        print(function_response)
        #存储图片,function_response是PIl.Image对象
        
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
        main_history.append({"role": "assistant", "content": second_response_message})
    else:
        main_history.append({"role": "assistant", "content": response_message["content"]})

    #解析输出，根据AI的指令调用函数

    output_text = []

    for i in range(0,len(main_history),2):
        output_text.append((main_history[i]['content'],main_history[i+1]['content']))

    return output_text,"2d_pic.png"



if __name__ == "__main__":
    with gr.Blocks() as demo:
        gr.Markdown("## MolecularGPT Chatbot")
        with gr.Row():
            with gr.Column():
                chat = gr.Chatbot().style(height=750)
                input_text = gr.Textbox(label='用户输入')
                button = gr.Button(label="提交")
            with gr.Column():
                image = gr.Image(interactive = False,value="2D_pic.png")
        button.click(chatbot,inputs= [input_text],outputs=[chat, image])
        input_text.submit(chatbot,inputs= [input_text],outputs=[chat, image])
    demo.launch()

