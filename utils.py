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
from matplotlib import colors
from rdkit.Chem.Draw import MolToImage
from torch_geometric.data import Data
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops
from rdkit.Chem import Descriptors
from torch_geometric.data import Data
import torch_geometric
from torch_geometric.data import DataLoader
from torch_geometric.nn import GCNConv

# import packages
# general tools
import numpy as np
# RDkit
from rdkit import Chem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
# Pytorch and Pytorch Geometric
import torch
from torch_geometric.data import Data
from torch.utils.data import DataLoader

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

# this function is used to save the image of the SMILES for vis purposes
def get_image(mol, atomset, name):
	hcolor = colors.to_rgb('green')
	if atomset is not None:
		#highlight the atoms set while drawing the whole molecule.
		img = MolToImage(mol, size=(600, 600),fitImage=True, highlightAtoms=atomset,highlightColor=hcolor)
	else:
		img = MolToImage(mol, size=(400, 400),fitImage=True)
	
	img = img.save("./media/" + name + ".jpg") 
	return img

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
        print ("adding node", i)
        print ("to drug", drug)
        G.add_node(i, attributes=drug)
    
    # Add edges to the graph
    for i, drug1 in enumerate(drug_data):
        print ("adding edge", i)
        print ("to drug", drug1)
        for j, drug2 in enumerate(drug_data):
            print ("adding edge",j)
            print ("to drug", drug2)
            if i != j:
                # Compute some similarity score between drug1 and drug2
                similarity = compute_similarity(drug1, drug2)
                threshold = 0.5
                if similarity > threshold:
                    # Add an edge between drug1 and drug2 if the similarity is above a certain threshold
                    G.add_edge(i, j, weight=similarity)
    
    return G

# convert the list of SMILES to a list of graphs
def smiles_to_graphs(smiles):
    graphs = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        adjacency = rdmolops.GetAdjacencyMatrix(mol)
        print (adjacency)
        x = torch.tensor([Descriptors.NumValenceElectrons(mol)], dtype=torch.float)
        edge_index = torch.tensor(adjacency.nonzero(), dtype=torch.long).t().contiguous()
        graphs.append(Data(x=x, edge_index=edge_index))
    return graphs

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
    
# code taken from https://www.blopig.com/blog/2022/02/how-to-turn-a-smiles-string-into-a-molecular-graph-for-pytorch-geometric/
    
def one_hot_encoding(x, permitted_list):
    """
    Maps input elements x which are not in the permitted list to the last element
    of the permitted list.
    """

    if x not in permitted_list:
        x = permitted_list[-1]

    binary_encoding = [int(boolean_value) for boolean_value in list(map(lambda s: x == s, permitted_list))]

    return binary_encoding

def get_atom_features(atom, use_chirality = True, hydrogens_implicit = True):
    """
    Takes an RDKit atom object as input and gives a 1d-numpy array of atom features as output.
    """

    # define list of permitted atoms
    
    permitted_list_of_atoms =  ['C','N','O','S','F','Si','P','Cl','Br','Mg','Na','Ca','Fe','As','Al','I', 'B','V','K','Tl','Yb','Sb','Sn','Ag','Pd','Co','Se','Ti','Zn', 'Li','Ge','Cu','Au','Ni','Cd','In','Mn','Zr','Cr','Pt','Hg','Pb','Unknown']
    
    if hydrogens_implicit == False:
        permitted_list_of_atoms = ['H'] + permitted_list_of_atoms
    
    # compute atom features
    
    atom_type_enc = one_hot_encoding(str(atom.GetSymbol()), permitted_list_of_atoms)
    
    n_heavy_neighbors_enc = one_hot_encoding(int(atom.GetDegree()), [0, 1, 2, 3, 4, "MoreThanFour"])
    
    formal_charge_enc = one_hot_encoding(int(atom.GetFormalCharge()), [-3, -2, -1, 0, 1, 2, 3, "Extreme"])
    
    hybridisation_type_enc = one_hot_encoding(str(atom.GetHybridization()), ["S", "SP", "SP2", "SP3", "SP3D", "SP3D2", "OTHER"])
    
    is_in_a_ring_enc = [int(atom.IsInRing())]
    
    is_aromatic_enc = [int(atom.GetIsAromatic())]
    
    atomic_mass_scaled = [float((atom.GetMass() - 10.812)/116.092)]
    
    vdw_radius_scaled = [float((Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) - 1.5)/0.6)]
    
    covalent_radius_scaled = [float((Chem.GetPeriodicTable().GetRcovalent(atom.GetAtomicNum()) - 0.64)/0.76)]

    atom_feature_vector = atom_type_enc + n_heavy_neighbors_enc + formal_charge_enc + hybridisation_type_enc + is_in_a_ring_enc + is_aromatic_enc + atomic_mass_scaled + vdw_radius_scaled + covalent_radius_scaled
                                    
    if use_chirality == True:
        chirality_type_enc = one_hot_encoding(str(atom.GetChiralTag()), ["CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW", "CHI_OTHER"])
        atom_feature_vector += chirality_type_enc
    
    if hydrogens_implicit == True:
        n_hydrogens_enc = one_hot_encoding(int(atom.GetTotalNumHs()), [0, 1, 2, 3, 4, "MoreThanFour"])
        atom_feature_vector += n_hydrogens_enc

    return np.array(atom_feature_vector)

def get_bond_features(bond, 
                      use_stereochemistry = True):
    """
    Takes an RDKit bond object as input and gives a 1d-numpy array of bond features as output.
    """

    permitted_list_of_bond_types = [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC]

    bond_type_enc = one_hot_encoding(bond.GetBondType(), permitted_list_of_bond_types)
    
    bond_is_conj_enc = [int(bond.GetIsConjugated())]
    
    bond_is_in_ring_enc = [int(bond.IsInRing())]
    
    bond_feature_vector = bond_type_enc + bond_is_conj_enc + bond_is_in_ring_enc
    
    if use_stereochemistry == True:
        stereo_type_enc = one_hot_encoding(str(bond.GetStereo()), ["STEREOZ", "STEREOE", "STEREOANY", "STEREONONE"])
        bond_feature_vector += stereo_type_enc

    return np.array(bond_feature_vector)

def create_pytorch_geometric_graph_data_list_from_smiles_and_labels(x_smiles, y):
    """
    Inputs:
    
    x_smiles = [smiles_1, smiles_2, ....] ... a list of SMILES strings
    y = [y_1, y_2, ...] ... a list of numerial labels for the SMILES strings (such as associated pKi values)
    
    Outputs:
    
    data_list = [G_1, G_2, ...] ... a list of torch_geometric.data.Data objects which represent labeled molecular graphs that can readily be used for machine learning
    
    """
    
    data_list = []
    
    for (smiles, y_val) in zip(x_smiles, y):
        
        # convert SMILES to RDKit mol object
        mol = Chem.MolFromSmiles(smiles)

        # get feature dimensions
        n_nodes = mol.GetNumAtoms()
        n_edges = 2*mol.GetNumBonds()
        unrelated_smiles = "O=O"
        unrelated_mol = Chem.MolFromSmiles(unrelated_smiles)
        n_node_features = len(get_atom_features(unrelated_mol.GetAtomWithIdx(0)))
        n_edge_features = len(get_bond_features(unrelated_mol.GetBondBetweenAtoms(0,1)))

        # construct node feature matrix X of shape (n_nodes, n_node_features)
        X = np.zeros((n_nodes, n_node_features))

        for atom in mol.GetAtoms():
            X[atom.GetIdx(), :] = get_atom_features(atom)
            
        X = torch.tensor(X, dtype = torch.float)
        
        # construct edge index array E of shape (2, n_edges)
        (rows, cols) = np.nonzero(rdmolops.GetAdjacencyMatrix(mol))
        torch_rows = torch.from_numpy(rows.astype(np.int64)).to(torch.long)
        torch_cols = torch.from_numpy(cols.astype(np.int64)).to(torch.long)
        E = torch.stack([torch_rows, torch_cols], dim = 0)
        
        # construct edge feature array EF of shape (n_edges, n_edge_features)
        EF = np.zeros((n_edges, n_edge_features))
        
        for (k, (i,j)) in enumerate(zip(rows, cols)):
            
            EF[k] = get_bond_features(mol.GetBondBetweenAtoms(int(i),int(j)))
        
        EF = torch.tensor(EF, dtype = torch.float)
        
        # construct label tensor
        y_tensor = torch.tensor(np.array([y_val]), dtype = torch.float)
        
        # construct Pytorch Geometric data object and append to data list
        data_list.append(Data(x = X, edge_index = E, edge_attr = EF, y = y_tensor))

    return data_list

# graph neural network model class
class GNN(torch.nn.Module):
    
    def __init__(self):
        super (GNN, self).__init__()
        self.conv1 = GCNConv(39, 64)

class GAE(nn.Module):
    def __init__(self, in_channels, out_channels):
        super(GAE, self).__init__()
        self.conv1 = GCNConv(in_channels, 2*out_channels)
        self.conv2 = GCNConv(2*out_channels, out_channels)
        self.conv3 = GCNConv(out_channels, in_channels)
        
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        
        x = self.conv1(x, edge_index)
        x = torch.relu(x)
        x = self.conv2(x, edge_index)
        x = torch.relu(x)
        x = self.conv3(x, edge_index)
        return x

      