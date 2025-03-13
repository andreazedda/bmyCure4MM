import torch
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from torch_geometric.nn import GAE, GCNConv
from torch_geometric.data import DataLoader
from colorama import Fore, init
import pickle

# Inizializza colorama per stampe colorate
init(autoreset=True)

# Definizione del modello Graph AutoEncoder
class GCNEncoder(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super(GCNEncoder, self).__init__()
        self.conv1 = GCNConv(in_channels, 2 * out_channels)
        self.conv2 = GCNConv(2 * out_channels, out_channels)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        return self.conv2(x, edge_index)

# Funzione per caricare i dati
def load_data(file_path):
    try:
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        print(Fore.GREEN + f"Dati caricati con successo da {file_path}")
        return data
    except Exception as e:
        print(Fore.RED + f"Errore durante il caricamento dei dati: {e}")
        return None

# Funzione di addestramento del modello
def train_model(data_list, in_channels, out_channels, epochs=100):
    model = GAE(GCNEncoder(in_channels, out_channels))
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    print(Fore.CYAN + "Inizio addestramento del modello...")
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()
        z = model.encode(data_list.x, data_list.edge_index)
        loss = model.recon_loss(z, data_list.edge_index)
        loss.backward()
        optimizer.step()

        if epoch % 10 == 0:
            print(Fore.YELLOW + f"Epoch {epoch}/{epochs} - Loss: {loss.item()}")

    print(Fore.GREEN + "Addestramento completato.")
    return model, z

# Funzione per la visualizzazione dello spazio latente
def plot_latent_space(Z, labels=None, apply_pca=True, n_clusters=3):
    Z_np = Z.detach().numpy()
    if apply_pca:
        pca = PCA(n_components=2)
        Z_np = pca.fit_transform(Z_np)
        print(Fore.YELLOW + "PCA applicato con successo.")

    if n_clusters:
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        clusters = kmeans.fit_predict(Z_np)
        plt.scatter(Z_np[:, 0], Z_np[:, 1], c=clusters, cmap='viridis', alpha=0.7)
        print(Fore.CYAN + "Cluster visualizzati con successo.")
    else:
        plt.scatter(Z_np[:, 0], Z_np[:, 1], alpha=0.7)

    plt.xlabel("Dimensione Latente 1")
    plt.ylabel("Dimensione Latente 2")
    plt.title("Visualizzazione dello Spazio Latente")
    plt.colorbar()
    plt.show()

# Pipeline principale
def bmyCure4MM_pipeline(file_path, in_channels, out_channels, epochs=100):
    data = load_data(file_path)
    if data is None:
        return

    model, Z = train_model(data, in_channels, out_channels, epochs)
    plot_latent_space(Z)

# Esempio di utilizzo
if __name__ == "__main__":
    file_path = "path_to_your_data.pkl"
    bmyCure4MM_pipeline(file_path, in_channels=50, out_channels=16, epochs=50)