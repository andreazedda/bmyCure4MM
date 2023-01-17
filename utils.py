#!pip install torch
#!pip install networkx
#!pip install rdkit
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import networkx as nx
import numpy as np
import rdkit.Chem as Chem
# drug_data=pd.read_csv("BBBP.csv")

# # Load the drug dataset (assume the data is in a CSV file with a column 'smiles')
# def load_drug_data():
#     drug_data = []
#     with open('BBBP.csv', 'r') as f:
#         for line in f:
#             smiles = line.strip()
#             drug_data.append({
#                 'smiles': smiles
#             })
#     return drug_data

