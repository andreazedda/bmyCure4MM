import pandas as pd


# Read in the data
drug_data=pd.read_csv("./BBBP.csv")

#sum all values equal to 1 in p_np column of drug_data
total_unitary_pnp = sum(drug_data.p_np==1)
print ("Total number of unitary p_np values: ", total_unitary_pnp)

#sum all values equal to 0 in p_np column of drug_data
total_zero_pnp = sum(drug_data.p_np==0)
print ("Total number of zero p_np values: ", total_zero_pnp)

# Print the total number of p_np values
print ("Total number of p_np values: ", len(drug_data.p_np))

# Calculate the percentage of unitary p_np values
print ("Percentage of unitary p_np values: ", round(total_unitary_pnp*100/len(drug_data.p_np),2), "%")

#calculate the percentage of zero p_np values
print ("Percentage of zero p_np values: ", round (total_zero_pnp*100/len(drug_data.p_np),2), "%")
from rdkit import Chem
from rdkit import RDLogger
from utils import *
from matplotlib import colors
from rdkit.Chem.Draw import MolToImage

# check the data frame
print(drug_data.head())


#extract smiles from the data frame and check the first 10 smiles
smiles = drug_data['smiles']
smiles.head()

#convert smiles to mols and disable warnings
RDLogger.DisableLog('rdApp.*')
mols = [Chem.MolFromSmiles(s) for s in smiles]
mols
#print (type(mols[0]))

#resizing the drug_data to
end_of_array = len (drug_data)
smiles=smiles[:end_of_array]
#print(smiles)

for i in range(0, len(smiles)):
    if smiles[i] == None:
        print(i)
# Create a graph representation of the first molecule
for i in range(0, end_of_array):
    name = drug_data['name'][i]
    if mols[i] != None:
        try: img = get_image(mols[i],None , name)
        except : print("un'immagine non è stata salvata")
from torch_geometric.data import DataLoader
import networkx as nx
from torch_geometric.nn import GAE
from torch_geometric.utils import train_test_split_edges
from utils import *


# Load the smiles and create the graph representation
smiles = drug_data['smiles'].to_list()
# Resize the array
smiles = smiles[:end_of_array]
res = [i for i in range(len(smiles)) if smiles[i] == '']
print (res)
# Load the labels we don't want to use them for training the autoencoder
labels = drug_data['p_np'].to_list()
# Resize the array
labels = labels[:end_of_array]
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
        # highlight the atoms set while drawing the whole molecule.
        img = MolToImage(mol, size=(600, 600), fitImage=True, highlightAtoms=atomset, highlightColor=hcolor)
    else:
        img = MolToImage(mol, size=(400, 400), fitImage=True)

    img = img.save("./media/" + name + ".jpg")
    return img




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


def get_atom_features(atom, use_chirality=True, hydrogens_implicit=True):
    """
    Takes an RDKit atom object as input and gives a 1d-numpy array of atom features as output.
    """

    # define list of permitted atoms

    permitted_list_of_atoms = ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe', 'As', 'Al', 'I',
                               'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'Li', 'Ge',
                               'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt', 'Hg', 'Pb', 'Unknown']

    if hydrogens_implicit == False:
        permitted_list_of_atoms = ['H'] + permitted_list_of_atoms

    # compute atom features

    atom_type_enc = one_hot_encoding(str(atom.GetSymbol()), permitted_list_of_atoms)

    n_heavy_neighbors_enc = one_hot_encoding(int(atom.GetDegree()), [0, 1, 2, 3, 4, "MoreThanFour"])

    formal_charge_enc = one_hot_encoding(int(atom.GetFormalCharge()), [-3, -2, -1, 0, 1, 2, 3, "Extreme"])

    hybridisation_type_enc = one_hot_encoding(str(atom.GetHybridization()),
                                              ["S", "SP", "SP2", "SP3", "SP3D", "SP3D2", "OTHER"])

    is_in_a_ring_enc = [int(atom.IsInRing())]

    is_aromatic_enc = [int(atom.GetIsAromatic())]

    atomic_mass_scaled = [float((atom.GetMass() - 10.812) / 116.092)]

    vdw_radius_scaled = [float((Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) - 1.5) / 0.6)]

    covalent_radius_scaled = [float((Chem.GetPeriodicTable().GetRcovalent(atom.GetAtomicNum()) - 0.64) / 0.76)]

    atom_feature_vector = atom_type_enc + n_heavy_neighbors_enc + formal_charge_enc + hybridisation_type_enc + is_in_a_ring_enc + is_aromatic_enc + atomic_mass_scaled + vdw_radius_scaled + covalent_radius_scaled

    if use_chirality == True:
        chirality_type_enc = one_hot_encoding(str(atom.GetChiralTag()),
                                              ["CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW",
                                               "CHI_OTHER"])
        atom_feature_vector += chirality_type_enc

    if hydrogens_implicit == True:
        n_hydrogens_enc = one_hot_encoding(int(atom.GetTotalNumHs()), [0, 1, 2, 3, 4, "MoreThanFour"])
        atom_feature_vector += n_hydrogens_enc

    return np.array(atom_feature_vector)


def get_bond_features(bond,
                      use_stereochemistry=True):
    """
    Takes an RDKit bond object as input and gives a 1d-numpy array of bond features as output.
    """

    permitted_list_of_bond_types = [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE,
                                    Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC]

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
        if mol is None:
            continue
        # get feature dimensions
        n_nodes = mol.GetNumAtoms()
        n_edges = 2 * mol.GetNumBonds()
        unrelated_smiles = "O=O"
        unrelated_mol = Chem.MolFromSmiles(unrelated_smiles)
        n_node_features = len(get_atom_features(unrelated_mol.GetAtomWithIdx(0)))
        n_edge_features = len(get_bond_features(unrelated_mol.GetBondBetweenAtoms(0, 1)))

        # construct node feature matrix X of shape (n_nodes, n_node_features)
        X = np.zeros((n_nodes, n_node_features))

        for atom in mol.GetAtoms():
            X[atom.GetIdx(), :] = get_atom_features(atom)

        X = torch.tensor(X, dtype=torch.float)

        # construct edge index array E of shape (2, n_edges)
        (rows, cols) = np.nonzero(rdmolops.GetAdjacencyMatrix(mol))
        torch_rows = torch.from_numpy(rows.astype(np.int64)).to(torch.long)
        torch_cols = torch.from_numpy(cols.astype(np.int64)).to(torch.long)
        E = torch.stack([torch_rows, torch_cols], dim=0)

        # construct edge feature array EF of shape (n_edges, n_edge_features)
        EF = np.zeros((n_edges, n_edge_features))

        for (k, (i, j)) in enumerate(zip(rows, cols)):
            EF[k] = get_bond_features(mol.GetBondBetweenAtoms(int(i), int(j)))

        EF = torch.tensor(EF, dtype=torch.float)

        # construct label tensor
        y_tensor = torch.tensor(np.array([y_val]), dtype=torch.float)

        # construct Pytorch Geometric data object and append to data list
        data_list.append(Data(x=X, edge_index=E, edge_attr=EF, y=y_tensor))

    return data_list


# graph neural network model class
class GNN(torch.nn.Module):

    def __init__(self):
        super(GNN, self).__init__()
        self.conv1 = GCNConv(39, 64)


class GAE(nn.Module):
    def __init__(self, in_channels, out_channels):
        super(GAE, self).__init__()
        self.conv1 = GCNConv(in_channels, 2 * out_channels)
        self.conv2 = GCNConv(2 * out_channels, out_channels)
        self.conv3 = GCNConv(out_channels, in_channels)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index

        x = self.conv1(x, edge_index)
        x = torch.relu(x)
        x = self.conv2(x, edge_index)
        x = torch.relu(x)
        x = self.conv3(x, edge_index)
        return x

data_list = create_pytorch_geometric_graph_data_list_from_smiles_and_labels(smiles, labels)
print ("First object in the data_list: " + str(data_list[0].x.shape))
import random
random.shuffle(data_list)
print ("Number of graphs in the data_list: " + str(len(data_list)))
print ("Number of graphs in the training set should be: " + str(int(len(data_list)*0.8)))
print ("Number of graphs in the test set should be: " + str(int(len(data_list)*0.2)))
train_len = int(len(data_list)*0.8)
test_len = int(len(data_list)*0.2)
train = data_list[:train_len] #train set
test = data_list[train_len:] #val set
print ("Number of graphs in the training set: " + str(len(train)))
print ("Number of graphs in the test set: " + str(len(test)))
from torch_geometric.loader import DataLoader
# imposing batch size
batch_size = len(train)
dataloader = DataLoader(dataset=train, batch_size=batch_size)


class GCNEncoder(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super(GCNEncoder, self).__init__()
        self.conv1 = GCNConv(in_channels, 2 * out_channels, cached=True)  # cached only for transductive learning
        self.conv2 = GCNConv(2 * out_channels, out_channels, cached=True)  # cached only for transductive learning

    # the forward function defines the computation performed at every call
    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        return self.conv2(x, edge_index)
from torch_geometric.nn import GCNConv, GAE
from torch_geometric.nn import GCNConv, GAE, GATConv, Linear, to_hetero


# parameters for the model
out_channels = 2
print ("The dimension of the output from the encoder part: " + str(out_channels))
num_features = 79
print ("Num of the features: " + str(num_features))

# Create the model
model = GAE(GCNEncoder(num_features, out_channels))
from torch_geometric.utils import scatter
from torch_geometric.loader import DataLoader
from torch.utils.tensorboard import SummaryWriter

epochs = 1000
print("Num of the epochs: " + str(epochs))

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

model = model.to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
graph_batch = next(iter(dataloader))

for epoch in range(1, epochs + 1):

    # Zero the gradients
    optimizer.zero_grad()

    graph_batch = graph_batch.to(device)

    batch_size = graph_batch.batch.max().item() + 1
    print("This is the batch size: " + str(batch_size))
    print(graph_batch)
    out = model.encode(graph_batch.x, graph_batch.edge_index)
    loss = model.recon_loss(out, graph_batch.edge_index)

    # Backward pass and optimization
    loss.backward()
    optimizer.step()

    # Logging
    if (i + 1) % 10 == 0:
        print("Epoch [{}/{}], Step [{}/{}], Loss: {:.4f}"
              .format(epoch + 1, 100, i + 1, len(dataloader), loss.item()))
Z= model.encode(graph_batch.x, graph_batch.edge_index)
print(Z)

import pandas as pd
import torch

# Load the BBBP dataset
drug_data = pd.read_csv("./BBBP.csv")

# Get the molecule IDs from the graph batch
graph_batch_molecule_ids = graph_batch.batch.unique()

# Filter the drug data to include only rows with molecule IDs in the graph batch
drug_data = drug_data[drug_data["num"].isin(graph_batch_molecule_ids)]

# Create a dictionary mapping each molecule ID to its label value
label_dict = dict(zip(drug_data["num"], drug_data["p_np"]))
from torch_geometric import utils
# Reindex the graph batch to make the batch indices zero-indexed and contiguous
graph_batch = utils.to_undirected(graph_batch)
graph_batch = utils.convert_node_labels_to_integers(graph_batch, node_attr="x")

# Convert the node features to a tensor
x = graph_batch.x.float()  # Convert to float tensor
edge_index = graph_batch.edge_index

# Create a list to store the labels for each atom
atom_labels = []

# Loop over each molecule in the graph batch
for num in graph_batch_molecule_ids:
    # Get the label value of the molecule
    molecule_label = label_dict[num]

    # Get the number of atoms in the molecule
    num_atoms = graph_batch.num_nodes[graph_batch.batch == num]

    # Assign the molecule's label to each atom
    atom_labels += [molecule_label] * num_atoms

# Load the embedding vectors from the numpy array
Z = model.encode(x, edge_index).detach().numpy()

# Convert the embedding vectors and labels to a pandas dataframe
z = pd.DataFrame(Z)
z["label"] = atom_labels

# Save the embedding vectors with the labels to a CSV file
z.to_csv("BBBP_embedding.csv", index=False)
