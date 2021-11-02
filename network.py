import numpy as np
import pandas as pd
import requests
import re
import sys
import matplotlib.pyplot as plt
import PIL
from PIL import Image
import io
import urllib.request
import itertools
import networkx as nx

# to display dataframe at full length
# will need to adjust width and number of max columns accordingly
desired_width=320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns',20)

def krest(keggid):
    """
    Download info from KEGG Rest API Module, process the info and create adjacent matrix from the info
    Args:
        keggid: Kegg ID of a Module in KEGG
    Returns:
        A_df: dataframe of the adjacent matrix
        G: networkx object of the adjacent matrix
    """
    # downloading the info from KEGG Rest API
    r = requests.get("http://rest.kegg.jp/link/rn/" + keggid)
    r2 = " ".join(r.text.split('\n'))
    r3 = re.split(r'(\t+)', r2)

    # putting the reactions in a list
    rl = []
    for line in r3:
        if line.startswith('rn'):
            rl.append(line)

    # unlisting list of reactions for later
    unlist_rl = " ".join(rl)

    # compiling reaction regex pattern
    rp = re.compile(r"\brn:\w{6}")

    # putting all reactions into one list
    result = rp.findall(unlist_rl)

    # grabbing chemical equation from reaction list and putting them into equation list (eql)
    eql = []
    for i in range(len(result)):
        url = "http://rest.kegg.jp/get/" + result[i]
        with requests.get(url) as file:
            lines = file.text.split('\n')
        for line in lines:
            if line.startswith('EQUATION '):
                eql.append(line)

    # splitting eql to substrate and product
    subs_l = []
    prod_l = []
    for i in range(len(eql)):
        subs = (eql[i].split('<='))[0]
        subs_l.append(subs)
        prod = (eql[i].split('<='))[1]
        prod_l.append(prod)

   # picking the compound id
    # setting a regex pattern for compound id
    patt = re.compile(r"\b\w{6}\b")
    subs_cpd = []
    for i in range(len(subs_l)):
        cpd1 = patt.findall(subs_l[i])
        subs_cpd.append(cpd1)

    prod_cpd = []
    for i in range(len(prod_l)):
        cpd2 = patt.findall(prod_l[i])
        prod_cpd.append(cpd2)

    # merging subs and prod compounds for each reaction
    rnl = []
    for i in range(len(eql)):
        rn = subs_cpd[i] + prod_cpd[i]
        rnl.append(rn)

    # creating enzyme unit matrix - if necessary
    eu = ['e{}'.format(j) for j in range(len(rnl))]

    # creating unique compound list for matrix construction
    unl = []
    for i in range(len(rnl)):
        unl.extend(rnl[i])

    un_set = set(unl)
    un = sorted(list(un_set))

    # combining uniq compound list and enzyme unit list to create matrix labels
    uneu = un + eu

    # creating adjacent matrix (adma) filled with zeroes
    mat_dim = len(un) + len(eu)
    A = np.zeros((mat_dim, mat_dim))
    # attaching labels to adma
    A_df = pd.DataFrame(A, index=uneu, columns=uneu)

    # assigning value of 1 (edge) to adma
    for i in range(len(rnl)):
        for j in range(len(rnl[i])):
            a = rnl[i][j]
            b = eu[i]
            A_df.loc[a, b] = 1
            A_df.loc[b, a] = 1

    # convert adma to networkx object
    G = nx.from_pandas_adjacency(A_df)
    return (A_df, G)


# function to determine nodes degree and centrality
def nodes_degree(df):
    G = nx.from_pandas_adjacency(df)
    d1 = dict(G.degree())
    centrality_degree = nx.degree_centrality(G)
    d2 = centrality_degree
    ds = [d1, d2]
    d = {}
    for k in d1.keys():
        d[k] = tuple(d[k] for d in ds)

    nodes_degree_df = pd.DataFrame.from_dict(d, orient='index', columns=['node degree', 'node centrality'])
    return nodes_degree_df


# function to work out the summary statistics of the nodes degree and/or centrality
def nodes_deg_summ(G):
    total_nodes = len(G)
    degree_values = dict(G.degree()).values()
    ave_degree = sum(d_values)/total_nodes
    sd_degree = np.std(d_values)
    med_degree = np.median(d_values)
    range_degree = np.max(d_values, axis=0) - np.min(d_values, axis=0)
    q3, q1 = np.percentile(d_values, [75 ,25])
    iqr_degree = q3 - q1
    cent_degree = nx.degree_centrality(G)
    data = {'average degree': [ave_degree],
            'std dev degree': [sd_degree],
            'median degree': [med_degree],
            'range degree': [range_degree],
            'IQR degree': [iqr_degree]}
    summary_df = pd.DataFrame(data)
    return summary_df

