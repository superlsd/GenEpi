import numpy as np
import networkx as nx
import itertools as it
import random as rd
import pickle as pk
import pandas as pd
import netmedpy

#Here, we define the PPI
ppi = pd.read_csv("input/autocore_ppi_symbol_lcc.csv",delimiter= ',',
           skipinitialspace=True)
G_ppi = nx.from_pandas_edgelist(ppi, 'symbol1', 'symbol2')

#Let's remove eventual nans and convert each string to index for speeding up the proces
node_to_index = {}
node_to_remove = []
i = 0
for n in G_ppi.nodes():
    if str(n)!='nan':
        node_to_index[n] = i
        i = i+1
    else:
        node_to_remove.append(n)

G_ppi.remove_nodes_from(node_to_remove)
G_ppi_relabel = nx.relabel_nodes(G_ppi, node_to_index)

#Let's import the pre-calculated spl dictionary
with open('input/ppi_spl.pickle', 'rb') as handle:
    spl = pk.load(handle)

#Let's convert the distance metrics into indexes
spl_index = {}
for node_pair,val in spl.items():
    if str(node_pair[0])!='nan' and str(node_pair[1])!='nan':
        spl_index[(node_to_index[node_pair[0]],node_to_index[node_pair[1]])] = val

#Let's import the diseases
with open('output/allsources_icd_gene_dict_unified.pickle', 'rb') as handle:
    icd_genes_dict = pk.load(handle)

#Let's convert the gene name lists into indexes
icd_genes_ppi_dict={}
for icd,genelist in icd_genes_dict.items():
    intersected_geneset = set(genelist)&set(G_ppi.nodes())
    intersected_geneset_index = [node_to_index[gene] for gene in intersected_geneset]
    if len(intersected_geneset_index)>1:
        icd_genes_ppi_dict[icd]=intersected_geneset_index

#Let's import the trajectories
with open('output/female_trajectories_dict_filtered_unique_agglomerated_genes.pickle', 'rb') as handle:
    female_trajectories_dict_filtered_unique_agglomerated_genes = pk.load(handle)

#Let's convert the gene name lists into indexes
female_trajectories_dict_filtered_unique_agglomerated_genes_idx={}
for trj,genelist in female_trajectories_dict_filtered_unique_agglomerated_genes.items():
    intersected_geneset = set(genelist)&set(G_ppi.nodes())
    intersected_geneset_index = [node_to_index[gene] for gene in intersected_geneset]
    if len(intersected_geneset_index)>1:
        female_trajectories_dict_filtered_unique_agglomerated_genes_idx[trj]=intersected_geneset_index

with open('output/male_trajectories_dict_filtered_unique_agglomerated_genes.pickle', 'rb') as handle:
    male_trajectories_dict_filtered_unique_agglomerated_genes = pk.load(handle)

#Let's convert the gene name lists into indexes
male_trajectories_dict_filtered_unique_agglomerated_genes_idx={}
for trj,genelist in male_trajectories_dict_filtered_unique_agglomerated_genes.items():
    intersected_geneset = set(genelist)&set(G_ppi.nodes())
    intersected_geneset_index = [node_to_index[gene] for gene in intersected_geneset]
    if len(intersected_geneset_index)>1:
        male_trajectories_dict_filtered_unique_agglomerated_genes_idx[trj]=intersected_geneset_index

#Let's import the drugs
with open('input/fda_approved_drug_gene_name_drugbank_2023.pickle', 'rb') as handle:
    fda_approved_drgs_dict = pk.load(handle)

fda_approved_drgs_ppi_dict = {}
for fda_drug,geneset in fda_approved_drgs_dict.items():
    intersected_geneset = set(geneset)&set(G_ppi.nodes())
    intersected_geneset_index = [node_to_index[gene] for gene in intersected_geneset]
    if len(intersected_geneset_index)>0:
        fda_approved_drgs_ppi_dict[fda_drug] = intersected_geneset_index

icd_fda_proximity_dict = {}
for icd,geneset1 in icd_genes_ppi_dict.items():
    for fda_drug,geneset2 in fda_approved_drgs_ppi_dict.items():
        prx_dict = netmedpy.proximity(G_ppi_relabel, geneset1,geneset2, spl_index,
                                          null_model="degree_match",n_iter=1000,
                                          symmetric=False)
        icd_fda_proximity_dict[icd,fda_drug] = prx_dict

with open('output/icd_fda_proximity_dict.pickle', 'wb') as handle:
    pk.dump(icd_fda_proximity_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

male_trj_fda_proximity_dict={}
for trj,geneset1 in male_trajectories_dict_filtered_unique_agglomerated_genes.items():
    for fda_drug,geneset2 in fda_approved_drgs_ppi_dict.items():
        prx_dict = netmedpy.proximity(G_ppi_relabel, geneset1,geneset2, spl_index,
                                          null_model="degree_match",n_iter=1000,
                                          symmetric=False)
        male_trj_fda_proximity_dict[trj,fda_drug] = prx_dict

with open('output/male_trj_fda_proximity_dict.pickle', 'wb') as handle:
    pk.dump(male_trj_fda_proximity_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

female_trj_fda_proximity_dict={}
for trj,geneset1 in female_trajectories_dict_filtered_unique_agglomerated_genes.items():
    for fda_drug,geneset2 in fda_approved_drgs_ppi_dict.items():
        prx_dict = netmedpy.proximity(G_ppi_relabel, geneset1,geneset2, spl_index,
                                          null_model="degree_match",n_iter=1000,
                                          symmetric=False)
        female_trj_fda_proximity_dict[trj,fda_drug] = prx_dict

with open('output/female_trj_fda_proximity_dict.pickle', 'wb') as handle:
    pk.dump(female_trj_fda_proximity_dict, handle, protocol=pk.HIGHEST_PROTOCOL)
