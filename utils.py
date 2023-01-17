import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import networkx as nx
import numpy as np
import rdkit.Chem as Chem
import networkx as nx
from rdkit import Chem
from rdkit import DataStructs
import torch.nn as nn

# Load the drug dataset (assume the data is in a CSV file with a column 'smiles')
def load_drug_data():
    drug_data = []
    with open('BBBP.csv', 'r') as f:
        for line in f:
            smiles = line.strip()
            drug_data.append({
                'smiles': smiles
            })
    return drug_data


#create a graph from a molecule
def mol_to_graph(mol):
    G = nx.Graph()
    G.add_nodes_from(range(mol.GetNumAtoms()))
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        G.add_edge(begin_atom, end_atom)
    return G

# Compute the similarity between two drugs based on their SMILES strings
# The similarity is computed using the Tanimoto coefficient
def compute_similarity(drug1, drug2):
    mol1 = Chem.MolFromSmiles(drug1)
    mol2 = Chem.MolFromSmiles(drug2)
    if mol1 is None or mol2 is None:
        return 0
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)
    return DataStructs.DiceSimilarity(fp1, fp2)


# Create a graph representation of the drug dataset
def create_graph(drug_data):
    # Create an empty graph
    G = nx.Graph()
    
    # Add nodes to the graph
    for i, drug in enumerate(drug_data):
        #print ("adding node", i)
        #print ("to drug", drug)
        G.add_node(i, attributes=drug)
    
    # Add edges to the graph
    for i, drug1 in enumerate(drug_data):
        #print ("adding edge", i)
        #print ("to drug", drug1)
        for j, drug2 in enumerate(drug_data):
            #print ("adding edge",j)
            #print ("to drug", drug2)
            if i != j:
                # Compute some similarity score between drug1 and drug2
                similarity = compute_similarity(drug1, drug2)
                threshold = 0.5
                if similarity > threshold:
                    # Add an edge between drug1 and drug2 if the similarity is above a certain threshold
                    G.add_edge(i, j, weight=similarity)
    
    return G

# Define the graph autoencoder model
class GraphAutoencoder(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(GraphAutoencoder, self).__init__()
        self.encoder = nn.Linear(input_dim, hidden_dim)
        self.decoder = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        x = self.encoder(x)
        x = self.decoder(x)
        return x