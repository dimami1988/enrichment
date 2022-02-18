import pandas as pd
import networkx as nx
from constants import CEND, CRED, EXPRESSION, INTERACTION, PATHWAY
import math
import numpy


def expressionLoader(expressionData):

    '''This function loads data from analyses of about 120 healthy and 20 cancer patients and sorts out most up or down regulated
    genes'''

    df = pd.read_csv(expressionData, sep="\t")
    downregulated = df[df.logFC < -0.5]
    upregulated = df[df.logFC > .5]
    result = pd.concat([upregulated, downregulated])
    result = result.dropna()
    print('\x1b[6;30;42m' + 'Data from research is loaded. Number of selected down/up regulated genes:' + '\x1b[0m')
    CRED = '\033[91m'
    CEND = '\033[0m'
    print(CRED + str(df.shape[0]) + CEND)
    print()
    return result


def convertLogFC(result):

    '''This function inverts logFC into binary 1 and -1'''

    col_one_list = result['Gene.symbol'].tolist()
    mask = (result['logFC'] < 0)
    result['logFC'] = result['logFC'].mask(mask, -1)
    mask = (result['logFC'] > 0)
    result['logFC'] = result['logFC'].mask(mask, 1)
    result = result.rename(columns={'Gene.symbol': 'symbol'})
    return result



def interactionLoader(interaction):

    '''This function loads a mapping file where there genes with up and down regulation
    that can be inverted into +1 and -1'''

    gene_interaction_map = pd.read_table(interaction)
    return gene_interaction_map


def graphCreate(pathway):
    '''This function loads  pathway and creates a graph out it'''
    df = pd.read_table(pathway, header=None)
    df = df.rename(columns={df.columns[0]: 'source'})
    df = df.rename(columns={df.columns[1]: 'weight'})
    df = df.rename(columns={df.columns[2]: 'target'})
    df = df[~df['target'].str.contains('CHEBI')]
    df = df[~df['source'].str.contains('CHEBI')]
    pathway = df
    print('\x1b[6;30;42m' + 'The length of pathway:' + '\x1b[0m')
    print(CRED + str(pathway.shape[0]) + CEND)
    print()
    #'''Here we reassign edges if they were found in mapping file'''
    G = nx.DiGraph()
    G = nx.from_pandas_edgelist(pathway, "source", "target", edge_attr=True)
    for node in G.nodes():
        G.nodes[node]['weight'] = 0
    return G


def addWeightsFromMap(G):

    '''This is gene mapping from somewhere so we can change weights in the path I found to inhibiton or activation'''

    for row in interactionLoader(INTERACTION).itertuples(index=True, name='Pandas'):
        if G.has_edge(row.source, row.target):
            G[row.source][row.target]['weight'] = row.relation
    return G

def edgeWeightAdd(G):

    '''This function iterates over Graph and assigns weights due to the type of interaction between genes'''

    for u, v, a in G.edges(data=True):
        if a['weight'] == 'activation':
            a['weight'] = 1
        if a['weight'] == 'inhibition':
            a['weight'] = -1
        else:
            a['weight'] = 0
    return G



def nodeWeightAdd(G, result):

    '''This function assigns weights to genes in graph which represents pathway, which is typical for the disease we research
     if they can be found in research analysis above.'''

    count = 0
    for row in result.index:
        if G.has_node(result.at[row, 'symbol']):
            count += 1
            G.nodes[result.at[row, 'symbol']]['weight'] = result.at[row, 'logFC']
    print(
        '\x1b[6;30;42m' + 'The number of genes from expression genes research analysis that appear in pathway:' + '\x1b[0m')
    print(CRED + str(count) + CEND)
    print()
    return G

def recalNodeWeights(G):
    '''This function calculates and updates weight of each node in graph by multiplying weight of edge by weight of node in HYP'''
    charge = 0
    for node in G.nodes():
        # print(node)
        leaves = list(G.adj[node])
        res = 0
        for el in leaves:
            # print(el)
            # print(G.nodes[node]['weight'])
            # print(G[node][el]["weight"])
            res += G.nodes[node]['weight'] * (G[node][el]["weight"])
        G.nodes[node]["weight"] = res
        # print(G.nodes[node]["weight"])
        # print(G.edges(node, data="weight"))
        # print()
    return G




def printing(G):
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



def concordance(G, result):
    '''This function calculates concordance'''
    p = 0.5
    problist = []
    hyplist = []
    for node in G.nodes():
        counter = 0
        li = 0
        ki = 0
        ni = 0

        leaves = list(G.adj[node])
        x = result['symbol'].isin(leaves)
        ni = sum(x)
        li = len(leaves) - sum(x)
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
    # print(probtable)
