from utils import *
import torch

num_features = 1
hidden_dim = 32
learning_rate = 0.01
num_epochs = 100
end_of_array = 10
drug_data=pd.read_csv("BBBP.csv")

# canonical training loop for a Pytorch Geometric GNN model gnn_model
gnn_model = GNN(num_features, hidden_dim).to(device)

smiles = drug_data['smiles'].to_list()
smiles = smiles[:end_of_array]
print (smiles)
labels = drug_data['p_np'].to_list()
labels = labels[:end_of_array]
print (labels)
data_list = create_pytorch_geometric_graph_data_list_from_smiles_and_labels(smiles, labels)
print (data_list[0])
# create dataloader for training
dataloader = DataLoader(dataset = data_list, batch_size = 2**7)

# define loss function
loss_function = nn.MSELoss()

# define optimiser
optimiser = torch.optim.Adam(gnn_model.parameters(), lr = 1e-3)

# loop over 10 training epochs
for epoch in range(10):

    # set model to training mode
    gnn_model.train()

    # loop over minibatches for training
    for (k, batch) in enumerate(dataloader):

        # compute current value of loss function via forward pass
        output = gnn_model(batch)
        loss_function_value = loss_function(output[:,0], torch.tensor(batch.y, dtype = torch.float32))

        # set past gradient to zero
        optimiser.zero_grad()

        # compute current gradient via backward pass
        loss_function_value.backward()

        # update model weights using gradient and optimisation method
        optimiser.step()