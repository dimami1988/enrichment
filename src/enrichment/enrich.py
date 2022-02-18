# https://github.com/cthoyt/cookiecutter-snekpack how to make folder structure
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-340 article for this analysis

import random
import pandas as pd
import networkx as nx

'''Test purpose graph. This function creates a graph out of random genes with randomly assigned edges 1 or -1 from the analysis above'''
#col_one_list = result['Gene.symbol'].tolist()
#G = nx.DiGraph()
#weights = [1, -1]
#for i in range(100):
#    nodes = random.sample(set(col_one_list), 2)
#    G.add_edge(nodes[0], nodes[1], weight=random.sample(set(weights), 1))

'''This function loads data from analyses of about 120 healthy and 20 cancer patients and sorts out most up or down regulated 
genes'''
#df=pd.read_csv(r'D:\Downloads\mechanism_enrichment\enrichment\data\GSE54129.top.table.tsv',sep="\t")
df=pd.read_csv(r'D:\Downloads\mechanism_enrichment\enrichment\data\GSE161533.top.table.tsv',sep="\t")
#print(df)
downregulated=df[df.logFC < -0.5]
upregulated=df[df.logFC > .5]
result = pd.concat([upregulated, downregulated])
result = result.dropna()
print('\x1b[6;30;42m' +   'Data from research is loaded. Number of selected down/up regulated genes:'+ '\x1b[0m')
CRED = '\033[91m'
CEND = '\033[0m'
print(CRED + str(df.shape[0])  + CEND)
print()

'''This function inverts logFC into binary 1 and -1'''
col_one_list = result['Gene.symbol'].tolist()
mask = (result['logFC'] <0)
result['logFC'] = result['logFC'].mask(mask, -1)
mask  = (result['logFC']  >0)
result['logFC'] = result['logFC'].mask (mask, 1)
result = result.rename(columns={'Gene.symbol': 'symbol'})



'''This function loads a mapping file where there genes with up and down regulation that can be inverted into +1 and -1'''
gene_interaction_map = pd.read_table( 'D:\Downloads\mechanism_enrichment\enrichment\data\gene_interaction_map.tsv')

'''This function loads signaling pathway and creates a graph out it'''
df=pd.read_table(r"D:\Downloads\VEGF_signaling_pathw.txt", header=None)
#df = pd.read_table(r"D:\Downloads\EGFR_Transactivation.txt", header=None)
#df = pd.read_table(r"D:\Downloads\Pentose_phosphate_pa.txt", header=None)
#df = pd.read_table(r"D:\Downloads\PI3K_AKT_activation.txt", header=None)
#df = pd.read_table(r"D:\Downloads\EGFR_downregulation.txt", header=None)
# Changing columns name with index number
df = df.rename(columns={df.columns[0]: 'source'})
df = df.rename(columns={df.columns[1]: 'weight'})
df = df.rename(columns={df.columns[2]: 'target'})
df = df[~df['target'].str.contains('CHEBI')]
df = df[~df['source'].str.contains('CHEBI')]
pathway=df
print('\x1b[6;30;42m' +   'The length of pathway:'+ '\x1b[0m')
print(CRED + str(pathway.shape[0])  + CEND)
print()

'''Here we reassign edges if they were found in mapping file'''
G = nx.DiGraph()
G = nx.from_pandas_edgelist(pathway, "source", "target", edge_attr=True)

for node in G.nodes():
    G.nodes[node]['weight']=0

'''This is gene mapping from somewhere so we can change weights in the path I found to inhibiton or activation'''
for row in gene_interaction_map.itertuples(index=True, name='Pandas'):
    if G.has_edge(row.source, row.target):
        G[row.source][row.target]['weight']  = row.relation

'''This function iterates over Graph and assigns weights due to the type of interaction between genes'''
for u, v, a in G.edges(data=True):
    if a['weight'] == 'activation':
        a['weight'] = 1
    if a['weight'] == 'inhibition':
        a['weight'] = -1
    else:
        a['weight'] = 0


#print(G.edges.data())
'''This function assigns weights to genes in graph which represents pathway, which is typical for the disease we research
 if they can be found in research analysis above.'''
count=0
for row in result.index:
    #print(result.at[row, 'symbol'] )
    if  G.has_node(result.at[row, 'symbol'] ):
        #print(result.at[row, 'symbol'] )
        #print(result.at[row, 'logFC'])
        count+=1
        G.nodes[result.at[row, 'symbol']]['weight']= result.at[row, 'logFC']
print('\x1b[6;30;42m' +   'The number of genes from expression genes research analysis that appear in pathway:'+ '\x1b[0m')
print(CRED + str(count)  + CEND)
print()


'''This function calculates and updates weight of each node in graph by multiplying weight of edge by weight of node in HYP'''
charge=0
for node in G.nodes():
    #print(node)
    leaves=list(G.adj[node])
    res=0
    for el in leaves:
        #print(el)
        #print(G.nodes[node]['weight'])
        #print(G[node][el]["weight"])
        res+= G.nodes[node]['weight'] * (G[node][el]["weight"])
    G.nodes[node]["weight"] = res
    #print(G.nodes[node]["weight"])
    #print(G.edges(node, data="weight"))
    #print()

'''This function shows sorted none zero values of graph'''
count = 0
genelist = []
weightlist = []
for node in G.nodes():
    if G.nodes[node]['weight'] != 0:
        # print(node)
        genelist.append(node)
        # print(G.nodes[node]['weight'])
        weightlist.append(G.nodes[node]['weight'])
        # print()
        count += 1

table = pd.DataFrame(
    {'gene': genelist,
     'weight': weightlist,
     })
table = table.sort_values(by=['weight'], ascending=False)
print(table)

'''This function calculates concordance'''
import math
import numpy

p = 0.5
# For some HYP x i
# ni # number of trials, n i , corresponds to the number of leaves that match gene expression data.
# ki # number of successful predictions k i is the number of leaves, mapped to state changes that are consistent with the HYP direction (correct).
# number of leaves of the same sign with node leaf

# li # Let l i be the number of downstream nodes for which the predicted direction cannot be determined (ambiguous) (not in gene expression data. ).

# probik= (math.factorial((ni-l))/(((math.factorial(k) * math.factorial(((n-l)-k)))   *  p**ki  *  (1−p)**(ni−ki−li)))


# subgraphs creating
problist = []
hyplist = []
for node in G.nodes():
    counter = 0
    li = 0
    ki = 0
    ni = 0

    leaves = list(G.adj[node])
    # x= result['symbol'].apply(lambda x: any([k in x for k in leaves]))
    x = result['symbol'].isin(leaves)
    ni = sum(x)
    li = len(leaves) - sum(x)
    # print(li)
    for el in leaves:
        signofleaf = numpy.sign(G.nodes[el]['weight'])
        signofnode = numpy.sign(G.nodes[node]['weight'])
        if signofleaf != 0 and signofleaf == signofnode:
            counter += 1
    ki = counter
    # print(ki)
    # print(ni-li-ki)
    power = ni - li - ki
    # print((ni-li))
    # (n-l)!/(k!((n-l)-k)!)
    x = (math.factorial(abs(ni - li))) / ((math.factorial(ki) * math.factorial(abs((abs(ni - li) - ki)))))

    probik = x * (p ** ki) * (1 - p) ** (power)
    problist.append(probik)
    hyplist.append(node)

probtable = pd.DataFrame(
    {'probability': problist,
     'HYP': hyplist,
     })
probtable = probtable.sort_values(by=['probability'], ascending=False)
pd.set_option('display.float_format', '{:.2f}'.format)
#print(probtable)