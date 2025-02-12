import numpy as np
import networkx as nx
import itertools as it
import pickle as pk
import pandas as pd

def aspdc(geneset_1,geneset_2,network):                #it returns the average among all the shortest path between two genesets
    spdc_dict = {}
    geneset_1 = set(set(geneset_1) & network.nodes())
    geneset_2 = set(set(geneset_2) & network.nodes())
    for i in geneset_1:
        for j in geneset_2:
            if (i,j) not in spdc_dict.keys() and (j,i) not in spdc_dict.keys():
                if i==j:
                    spdc_dict[(i,j)]=0
                else:
                    try:
                        try:
                            spdc_dict[(i,j)] = spl[i,j]
                        except:
                            spdc_dict[(i,j)] = spl[j,i]
                    except:
                        spdc_dict[(i,j)]=nx.shortest_path_length(network,i,j)
    aspdc_value = sum(list(spdc_dict.values()))/len(spdc_dict)
    return aspdc_value

def mean_shortest_distance_single(geneset_1,network):                #it returns the dab
    spdc_dict = {}
    all_distances = []
    geneset_1 = set(set(geneset_1) & network.nodes())
    for geneA in geneset_1:
        all_distances_A = []
        for geneB in geneset_1:
            if (geneA,geneB) not in spdc_dict.keys() and (geneB,geneA) not in spdc_dict.keys():
                if geneA!=geneB:     #We exclude those genes that are the same
                    try:
                        try:
                            sh_path = spl[geneA,geneB]
                            spdc_dict[(geneA,geneB)]= sh_path
                            all_distances_A.append(sh_path)
                        except:
                            sh_path = spl[geneB,geneA]
                            spdc_dict[(geneA,geneB)]= sh_path
                            all_distances_A.append(sh_path)
                    except:
                        sh_path = nx.shortest_path_length(network,geneA,geneB)
                        spdc_dict[(geneA,geneB)]= sh_path
                        all_distances_A.append(sh_path)
        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)
    mean_shortest_distance = np.mean(all_distances)
    return mean_shortest_distance


def mean_shortest_distance_pair(geneset_1,geneset_2,network):
    gene_set1 = set(geneset_1) & network.nodes()
    gene_set2 = set(geneset_2) & network.nodes()
    spdc_dict = {}
    all_distances = []
    for geneA in gene_set1:
        all_distances_A = []
        for geneB in gene_set2:
            if (geneA,geneB) not in spdc_dict.keys() and (geneB,geneA) not in spdc_dict.keys():
                if geneA==geneB:   #Between sets we keep also the overlapping genes
                    spdc_dict[(geneA,geneB)]=0
                    all_distances_A.append(0)
                else:
                    try:
                        try:
                            sh_path = spl[geneA,geneB]
                            spdc_dict[(geneA,geneB)]= sh_path
                            all_distances_A.append(sh_path)
                        except:
                            sh_path = spl[geneB,geneA]
                            spdc_dict[(geneA,geneB)]= sh_path
                            all_distances_A.append(sh_path)
                    except:
                        sh_path = nx.shortest_path_length(network,geneA,geneB)
                        spdc_dict[(geneA,geneB)]= sh_path
                        all_distances_A.append(sh_path)
        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)
    for geneA in gene_set2:
        all_distances_A = []
        for geneB in gene_set1:
            if (geneA,geneB) not in spdc_dict.keys() and (geneB,geneA) not in spdc_dict.keys():
                if geneA==geneB:
                    spdc_dict[(geneA,geneB)]=0
                    all_distances_A.append(0)
                else:
                    try:
                        try:
                            sh_path = spl[geneA,geneB]
                            spdc_dict[(geneA,geneB)]= sh_path
                            all_distances_A.append(sh_path)
                        except:
                            sh_path = spl[geneB,geneA]
                            spdc_dict[(geneA,geneB)]= sh_path
                            all_distances_A.append(sh_path)
                    except:
                        sh_path = nx.shortest_path_length(network,geneA,geneB)
                        spdc_dict[(geneA,geneB)]= sh_path
                        all_distances_A.append(sh_path)
        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)
    # calculate mean shortest distance
    mean_shortest_distance = np.mean(all_distances)
    return mean_shortest_distance

def separation(geneset_1,geneset_2,network):                #it returns the separation coefficient between two genesets
    d_ab = mean_shortest_distance_pair(geneset_1,geneset_2,network)
    d_aa = mean_shortest_distance_single(geneset_1,network)
    d_bb = mean_shortest_distance_single(geneset_2,network)
    s_ab = d_ab - (d_aa+d_bb)/2
    return s_ab

def overlap_jaccard(list1,list2):                          #it returns the jaccard index
    intersction_term= len(set(list1) & set(list2))
    denominator = len(set(list1).union(set(list2)))
    overlap_jaccard_coeff = intersction_term/denominator
    return overlap_jaccard_coeff

def overlap_coefficient(list1,list2):                          #it returns the jaccard index
    intersction_term= len(set(list1) & set(list2))
    denominator = min([len(set(list1)),len(set(list2))])
    overlap_coeff = intersction_term/denominator
    return overlap_coeff

def calculate_closest_distance(G_ppi_lcc,spl, nodes_from, nodes_to):
    values_outer = []
    for node_from in nodes_from:
        values = []
        for node_to in nodes_to:
            if node_from==node_to:
                val =0
            else:
                try:
                    try:
                        val = spl[node_from,node_to]
                    except:
                        val = spl[node_to,node_from]
                except:
                    val=len(nx.shortest_path(G_ppi_lcc,source=node_from, target=node_to))
            values.append(val)
        d = min(values)
        #print d,
        values_outer.append(d)
    d = np.mean(values_outer)
    #print d
    return d

#Here, we define the PPI
ppi = pd.read_csv("input/autocore_ppi_symbol_lcc.csv",delimiter= ',',
           skipinitialspace=True)
G_ppi = nx.from_pandas_edgelist(ppi, 'symbol1', 'symbol2')
print("The total PPI is read and it is big: %s" %(len(G_ppi)))

#Let's import the pre-calculated spl dictionary
with open('input/ppi_spl.pickle', 'rb') as handle:
    spl = pk.load(handle)

#Let's import the Elma diseases
with open('output/allsources_icd_gene_dict_unified.pickle', 'rb') as handle:
    icd_genes_dict = pk.load(handle)

icd_genes_ppi_dict={}
for icd,genelist in icd_genes_dict.items():
    newgenelist=set(genelist)&set(G_ppi.nodes())
    if len(newgenelist)>1:
        icd_genes_ppi_dict[icd]=newgenelist

pairwise_diseases= list(it.combinations(list(icd_genes_ppi_dict.keys()), 2))
icd_genes_ppi_overlap_coefficient_dict = {}
icd_genes_ppi_overlap_jaccard_dict = {}
icd_genes_ppi_aspdc_dict = {}
icd_genes_ppi_mean_shortest_distance_dict = {}
icd_genes_ppi_separation_dict = {}

dis_list = list(icd_genes_ppi_dict.keys())

for dis_pair in pairwise_diseases:
    geneset1=icd_genes_ppi_dict[dis_pair[0]]
    geneset2=icd_genes_ppi_dict[dis_pair[1]]
    icd_genes_ppi_overlap_coefficient_dict[dis_pair] = overlap_coefficient(geneset1,geneset2)
    icd_genes_ppi_overlap_jaccard_dict[dis_pair] = overlap_jaccard(geneset1,geneset2)
    icd_genes_ppi_aspdc_dict[dis_pair] = aspdc(geneset1,geneset2,G_ppi)
    icd_genes_ppi_mean_shortest_distance_dict[dis_pair] = mean_shortest_distance_pair(geneset1,geneset2,G_ppi)
    icd_genes_ppi_separation_dict[dis_pair] = separation(geneset1,geneset2,G_ppi)

with open('output/icd_genes_ppi_overlap_coefficient_dict.pickle', 'wb') as handle:
    pk.dump(icd_genes_ppi_overlap_coefficient_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

with open('output/icd_genes_ppi_overlap_jaccard_dict.pickle', 'wb') as handle:
    pk.dump(icd_genes_ppi_overlap_jaccard_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

with open('output/icd_genes_ppi_aspdc_dict.pickle', 'wb') as handle:
    pk.dump(icd_genes_ppi_aspdc_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

with open('output/icd_genes_ppi_mean_shortest_distance_dict.pickle', 'wb') as handle:
    pk.dump(icd_genes_ppi_mean_shortest_distance_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

with open('output/icd_genes_ppi_separation_dict.pickle', 'wb') as handle:
    pk.dump(icd_genes_ppi_separation_dict, handle, protocol=pk.HIGHEST_PROTOCOL)
