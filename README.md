# KEGG Network Analysis

Python utilities for retrieving reaction information from the KEGG REST API and representing metabolic modules as networks.

The functions in this repository were developed for exploratory pathway and network analysis in biomedical research.

## What it does

- retrieves reactions associated with a KEGG module
- extracts compounds from reaction equations
- constructs a compound–reaction adjacency matrix
- converts the matrix to a NetworkX graph
- calculates node degree and degree centrality
- summarises degree characteristics across the network
