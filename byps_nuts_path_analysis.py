import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import itertools
from tqdm import tqdm

#-------------------------------------------------------------------------
# Generating the bipartite graph.
#-------------------------------------------------------------------------
def genScopeGraph(scopeRxns, rxnMat, prodMat):
    # Calculating the scopeMets.
    scopeRxnVec = np.zeros(len(rxnMat))
    scopeRxnVec[scopeRxns] = 1
    scopeMetVec = (np.dot(np.transpose(rxnMat), scopeRxnVec) > 0) * 1
    scopeMetVec = (np.dot(np.transpose(prodMat), scopeRxnVec) + scopeMetVec > 0) * 1
    scopeMets = np.nonzero(scopeMetVec)[0]

    # Trying to write a bipartite graph of the expanded scope.
    scopeGraph = nx.DiGraph()
    scopeGraph.add_nodes_from([id_to_kegg[x] for x in scopeMets], bipartite = 0)
    scopeGraph.add_nodes_from([rxn_id_to_kegg[x] for x in scopeRxns], bipartite = 1)

    # Connect all the metabolites to reactions.
    for r in scopeRxns:
        scopeGraph.add_edges_from([(id_to_kegg[c], rxn_id_to_kegg[r])
                                    for c in np.nonzero(rxnMat[r])[0]])
        scopeGraph.add_edges_from([(rxn_id_to_kegg[r], id_to_kegg[p]) 
                                    for p in np.nonzero(prodMat[r])[0]])

    return scopeGraph


def giveScopeRxns(sources, thisGenotype):
    kos = [pos_to_gene_map[i] for i in np.where((thisGenotype > 0.5) * 1)[0]]
    refDF = pd.read_csv('KOREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])
    def give_number( string ):
        return int( string[-5:] )
    refDF['id'] = refDF['id'].apply( give_number )
    refDF['rid'] = refDF['rid'].apply( give_number )
    can = refDF.loc[refDF['id'].isin(kos),:]['rid'].values
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in can]
    can = list((set(can)) & set(rxns))
    avrxns = np.array([rxn_kegg_to_id[e] for e in can])
    scopeMets, avScopeRxn = giveScope(rxnMat[avrxns], prodMat[avrxns], seedVec, 
                                     sumRxnVec[avrxns])
    scopeRxns = avrxns[np.nonzero(avrxns * avScopeRxn)[0]]
    return scopeRxns


def giveAddAnc(anc, gainedSet):
    addAnc = np.copy(anc.genotype)
    for p in gainedSet:
        addAnc[gene_to_pos_map[p]] = 1

    return addAnc

def giveYield(path):
    mulFac = [1]
    energy_cpds = ['C00002']
    Energy = [kegg_to_id[i] for i in energy_cpds]
    rxns_s = [el for i, el in enumerate(path) if i % 2 == 1]
    rxns = np.array([rxn_kegg_to_id[i] for i in rxns_s])
    return int( sum( [ sum( stoich_matrix[ rxn ][ Energy ] * mulFac ) 
                for rxn in rxns.astype(int) ] ) )

# Figuring out the path lengths from HGT byproducts to core.
lengths_nuts, yields_nuts = [], []
lengths_byps, yields_byps = [], []
for thisNode in tqdm(nodes):
    for thisChild in thisNode.children:
        gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
        if gainGenes:
            des = giveAddAnc(thisNode, gainGenes)
            scopeRxns = giveScopeRxns(ubyps, des)
            scopeGraph = genScopeGraph(scopeRxns, rxnMat, prodMat)
            nscopeRxns = giveScopeRxns(unuts, des)
            nscopeGraph = genScopeGraph(nscopeRxns, rxnMat, prodMat)
            for tbyp in tqdm(ubyps):
                if tbyp in scopeGraph:
                    paths = nx.shortest_path_length(scopeGraph, tbyp)
                    if len(paths) > 1:
                        for cpd in paths:
                            if cpd in Core:
                                lengths_byps.append(paths[cpd])
                                yields_byps.append(giveYield(
                                                   nx.shortest_path(scopeGraph, tbyp, cpd)))
                                lengths_nuts.append([])
                                yields_nuts.append([])
                                for tnut in unuts:
                                    if tnut in nscopeGraph:
                                        npaths = nx.shortest_path_length(nscopeGraph, tnut)
                                        if len(npaths) > 1:
                                            if cpd in npaths:
                                                    lengths_nuts[-1].append(npaths[cpd])
                                                    yields_nuts[-1].append(giveYield(
                                                        nx.shortest_path(nscopeGraph, tnut, cpd)))

diff_lengths = []
for i, tl in enumerate(lengths_byps):
    if lengths_nuts[i]:
        diff_lengths.append(tl - np.median(lengths_nuts[i]))
    else:
        continue

sns.kdeplot(np.asarray(diff_lengths), bw=1)
plt.show()

diff_yields = []
for i, tl in enumerate(yields_byps):
    if yields_nuts[i]:
        diff_yields.append(tl - np.mean(yields_nuts[i]))
    else:
        continue

sns.kdeplot(np.asarray(diff_yields), bw=1)
plt.show()


lengths_nuts_u = unlistify(lengths_nuts)
yields_nuts_u = unlistify(yields_nuts)

sns.kdeplot(np.asarray(lengths_byps), bw=1, color='r')
sns.kdeplot(np.asarray(lengths_nuts_u), bw=1, color='grey')
plt.show()
