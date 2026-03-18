import networkx as nx
import numpy as np
 
def build_ppi_network():
    G=nx.Graph()
    proteins=['TP53','MDM2','BRCA1','ATM','CDK2','CCND1','RB1','E2F1','MYC','MAX','BCL2','BAX']
    interactions=[('TP53','MDM2'),('TP53','BRCA1'),('TP53','ATM'),('CDK2','CCND1'),
                   ('CDK2','RB1'),('RB1','E2F1'),('MYC','MAX'),('TP53','MYC'),('BCL2','BAX'),('MDM2','BCL2')]
    for p in proteins: G.add_node(p)
    G.add_edges_from(interactions)
    return G
 
def find_hub_proteins(G, top_k=3):
    deg=dict(G.degree()); return sorted(deg,key=deg.get,reverse=True)[:top_k]
 
def find_modules(G):
    return list(nx.connected_components(G))
  def betweenness_centrality(G):
    return nx.betweenness_centrality(G)
 
def essential_proteins(G, threshold=0.1):
    bc=betweenness_centrality(G)
    return [p for p,s in bc.items() if s>=threshold]
 
G=build_ppi_network()
print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
print(f"Hub proteins (top-3): {find_hub_proteins(G)}")
print(f"Connected modules: {len(find_modules(G))}")
print(f"Essential proteins (BC≥0.1): {essential_proteins(G)}")
try: print(f"Clustering coefficient: {nx.average_clustering(G):.3f}")
except: pass
