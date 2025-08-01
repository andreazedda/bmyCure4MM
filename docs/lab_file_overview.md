# Lab Folder Overview

This document provides brief explanations of the scripts and notebooks stored in the `lab` directory. These files are mainly prototypes and exploratory experiments used during research and testing of the project.

## Contents

### Notebooks

- **GAE_on_multiple_graphs_in _multiple_batch.ipynb** – Prototype of a graph
  autoencoder trained on batches of multiple molecules. The notebook loads the
  BBBP dataset and attempts to reconstruct each molecule from its graph
  representation.
- **GAE_on_multiple_graphs_in _single_batch.ipynb** – Similar experiment where
  all graphs are processed within a single batch.
- **GAE_on_multiple_graphs_in _single_batch copy.ipynb** – Duplicate of the
  single‑batch experiment kept for reference.
- **GAE_on_multiple_graphs_old.ipynb** – Early version of the multi‑graph
  autoencoder workflow.
- **GAE_on_single_graph.ipynb** – Exploratory analysis of the BBBP dataset with
  background information about the dataset and an attempt at building a graph
  autoencoder for a single molecule.
- **GDL_on_demographics.ipynb** – Loads a public hospital demographics CSV file
  and demonstrates basic graph deep learning operations on the data.
- **GDL_on_demographics_MM.ipynb** – Variation of the previous notebook using a
  multiple myeloma specific dataset.
- **building_graph.ipynb** – Example showing how to construct a heterogeneous
  graph from the MovieLens dataset and prepare train/validation/test splits with
  PyTorch Geometric.
- **converting_data_from_gudhi_to_PG.ipynb** – Small snippet illustrating how to
  translate a Gudhi simplicial complex into a `torch_geometric` graph object.
- **enzymes.ipynb** – Short explanation of PyTorch Geometric `DataLoader`
  behaviour, used as a quick reference.
- **graph_plotting.ipynb** – Demonstrates 3‑D plotting of a simple graph using
  NetworkX and Matplotlib.
- **graphs_classifier.ipynb** – Notebook exploring graph classification for the
  BBBP dataset.
- **load_csv_example.ipynb** – Despite the extension, this file contains pure
  Python code that downloads the MovieLens CSV files and converts them into
  `HeteroData` objects.
- **solubility_GDL.ipynb** – Prototype applying geometric deep learning methods
  to predict drug solubility.

### Python scripts

- **Autoencoder.py** – Loads the BBBP dataset, converts SMILES to RDKit
  molecules and prints basic statistics as a starting point for building a graph
  autoencoder.
- **latent_drug_map.py** – Pipeline for training a `GAE` model on prepared
  graphs and visualising the latent space with PCA and K‑means clustering.
- **pytorch_mlp_example.py** – Minimal implementation of a multilayer perceptron
  trained on the CIFAR‑10 image dataset.
- **training.py** – Example training loop for a `GNN` model using molecular
  graphs generated from the BBBP dataset.
- **utils.py** – Collection of helper functions for working with SMILES strings,
  converting them into graph objects and building Gudhi complex files.

