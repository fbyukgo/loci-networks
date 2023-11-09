#!/usr/bin/env python3
#packages
import sys
import os.path
import csv
from igraph import *
import time
import pandas as pd
from itertools import chain
from itertools import combinations
import argparse
#args, taking .clumped files

start1 = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--types', nargs='+')
parser.add_argument('--files', nargs='+')
parser.add_argument('--analysis_name')
parser.add_argument('--summary_dir')

analysis_list = parser.parse_args().types
f_list = parser.parse_args().files
analysis_name = parser.parse_args().analysis_name
summary_dir = parser.parse_args().summary_dir

file_dict = {}
for _ in range(len(analysis_list)):
    typ = analysis_list[_]
    f = f_list[_]
    file_dict[typ] = f

print(file_dict)
summary_full_path = os.path.join(summary_dir, "%s.txt" % analysis_name)


#functions
# find each analysis combination from given analysis list and set values
#for example for bottomline only analysis it sets:
#bottomline={bottomline:1, minp:0, naive:0, largest:0} in the output dictionary
def get_analysis_combination(analysis_list):
    combination_dict = {}
    analysis_dict = {}
    c = 1
    while c <= len(analysis_list):
        combination = list(combinations(analysis_list, c))
        combination_dict[c] = combination
        c += 1

    for comb  in combination_dict.values():
        for ele in comb:
            analysis_name = "_".join(ele)
            status_dict = {}
            for analysis in analysis_list:
                status_dict[analysis] = 0
                if analysis in ele:
                    status_dict[analysis] = 1
                analysis_dict[analysis_name] = status_dict
    return analysis_dict

#store index SNPs for given file, currently unused
def store_index_snps(f):
    index_snps = []

    with open(f) as file_in:
        for line in file_in.read().splitlines():
            fields = line.split()
            index_snp = fields[2]
            index_snps.append(index_snp)


    return list(set(index_snps))


#make graph for each analysis file and store them seperately for each chromosome in output dictionary(unioned_graphs)
#the function takes union of  them before storing in order merger clumps from different ancestries
def graph_maker(f, analysis):
    graphs = {}
    unioned_graphs = {}
    for chrm in range(1, 23):
        graphs[chrm] = []
        unioned_graphs[chrm] = []
    with open(f) as file_in:
        for line in file_in.read().splitlines():
            fields = line.split()
            chromosome = int(fields[0])
            index_snp = fields[2]
            other_snps = fields[11]

            if other_snps == "NONE":
                L_tmp=[]
                L_tmp.append(index_snp)
            else:
                L_tmp = [ snp.split("(")[0] for snp in other_snps.split(",")]
                L_tmp.append(index_snp)

            L=[]
            for snp in L_tmp:
                if snp not in L:
                    L.append(snp)
            G = Graph.Full(n=len(L))
            G.vs["name"] = L
            G.vs["group"] = analysis
            graphs[chromosome].append(G)


    for k, v in graphs.items():
        if v:
            U = union(v, byname="auto")
            unioned_graphs[k] = U

    #U = union(graphs, byname="auto")
    #U = disjoint_union(graphs)
    return unioned_graphs


# not in use but can be used in order to retrieve a cluster graph or number of connected components from a given graph
def get_connected_components(graph, r):
    clusters = graph.clusters( mode = "weak")
    ncomponents = len(clusters)

    if r == "number":
        return ncomponents
    elif r == "graph":
        return clusters


#get union across the analysis types
def get_union(graphs_dict):
    cluster_graphs = []

    for k in range(1, 23):
        graphs = []
        for analysis in graphs_dict.keys():
            graphs.append(graphs_dict[analysis][k])
        if graphs:
            graphs_no_empty = [ graph for graph in graphs if graph != []]
            print(f"for chr {k} there are {len(graphs)} graphs to combine")
            #graphs = minp_graph + bottomline_graph + largest_graph
            union_graph=union(graphs_no_empty, byname='auto')
            cluster_graph = get_connected_components(union_graph, "graph")
            cluster_graphs.append(cluster_graph)

    return cluster_graphs


def intersection(lst1, lst2):
        return  list(set(lst1) & set(lst2))

#find exclusive components based on the "group" vertex attribute
def find_exclusive_components(analysis_list, analysis_dict, cluster_graphs):

    count_dict = {}
    for k in analysis_dict.keys():
        count_dict[k] = 0

    for cluster_graph in cluster_graphs:
        for i in range(len(cluster_graph)):

            subgraph = cluster_graph.subgraph(i)
            subgraph_array = subgraph.vs["name"]

            attr_list = []
            for attr in subgraph.vs.attributes():
                analysis = subgraph.vs[attr]
                attr_list.append(analysis)

            attr_list_flat = list(chain(*attr_list))


            status_dict = {}
            #set all analysis intersections 0
            for analysis in analysis_list:
                status_dict[analysis] = 0
                if analysis in attr_list_flat:
                    status_dict[analysis] = 1

            for k, v in analysis_dict.items():
                if v == status_dict:
                    count_dict[k] += 1


            #print(f"Signal {i} subgraph:")
            #print(subgraph.vs["name"])

    str_list = []
    for k, v in count_dict.items():
        count_string = f"{k}: {v}"
        str_list.append(count_string)


    final_string = '\n'.join([*str_list])
    with open(summary_full_path, 'w') as f:
        print(final_string, file=f)

load = time.time()

print(f"all packages loaded and functions are defined in {load - start1} second")

#functions in use

#1 make graphs for both analysis sets

#making graphs for each analysis type and getting internal union
graphs_dict = {}
for name, filedir in file_dict.items():
    print(f"running graph make on {name}")
    graph = graph_maker(filedir, name)
    graphs_dict[name] = graph


#print("gettin index snps for minp")
#minp_index_snps = store_index_snps(minp_file)

graph_made = time.time()
print(f"all graphs are made in {graph_made - load} second")


#2 union  graphs and generate cluster graph

cluster_graphs = get_union(graphs_dict)
union_made = time.time()
print(f"the union graph is generated in {union_made - graph_made} seconds")

# generate all analysis combinations from given analysis list
analysis_dict = get_analysis_combination(analysis_list)

#3 find exclusive component(loci)
find_exclusive_components(analysis_list, analysis_dict,  cluster_graphs)

end1 = time.time()
print(f"exclusive components are found in {end1 - union_made} seconds")
print(f"whole script takes { end1 - start1 } second to run")



