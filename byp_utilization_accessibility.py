from ete3 import Tree
import numpy as np
import pandas as pd
from uniqify import uniqify
from unlistify import unlistify
from copy import copy
import pickle
import re
from load_kegg import *

# Getting the set of bacterial abbreviations for each organism in KEGG.
orgNames = []
with open('endo_removed_prok_abbr_kegg.txt', 'r') as f:
    for thisLine in f.readlines():
        orgNames.append( thisLine.strip() )


geneDict = {thisOrg : 
            set(np.genfromtxt('organism_kogenes/' + thisOrg + '.txt').astype(int))
            for thisOrg in orgNames}
gene_ids = sorted( uniqify( unlistify( list( geneDict.values() ) ) ) )

# Getting the tree with internal nodes from gainLoss' ancestral reconstruction.
# Internal nodes are labeled with 'Nx' where x is a number, the root being
# '[N1]' and then the numbers increase.
gl_tree = Tree( 'full_proks_gainLoss_results/TheTree.INodes.ph', format=1 )
ind_tree = Tree('ancs_fastanc_indscores.txt', format = 1)

# # To print the tree.
# print(gl_tree.get_ascii(show_internal=True))

# Getting the root of the tree.
root = gl_tree&'[N1]'
nodes = list( root.traverse() )

# Getting the original traits in a trait dict.
ind_degree = pickle.load(open('ind_degree_dict.dat', 'rb'))
kegg_to_tip_name_dict = pickle.load(open('kegg_to_tip_name_dict.dat', 'rb'))
tip_to_name_segata_dict = pickle.load(open('tip_to_name_segata_dict.dat', 'rb'))
for orgName in ind_degree:
    full_name = tip_to_name_segata_dict[kegg_to_tip_name_dict[orgName]]    
    (gl_tree&full_name).add_feature( 'ind_degree', len(ind_degree[ orgName ]) )

# int(float(find_between_r(nodes[1].name, '{', '}')))

# # Writing all internal ancestral states.
# def markChildren(gl_node, anc_node):
#     gl_node.add_feature('ind_degree', int(float(re.findall(r'\{(.*?)\}', anc_node.name)[0])))

#     for gl_child, anc_child in zip(gl_node.children, anc_node.children):
#         try:
#             gl_child.ind_degree
#         except:
#             markChildren(gl_child, anc_child)

# markChildren(root, ind_tree)

# Getting the original genotypes in a genotype dict.
# NOTE: genotype dicts are now probability vectors.
# for orgName in isIndDict:
#         (gl_tree&orgName).add_feature( 'genotype', 
#                                 np.array( list( map( int, genotype_dict[ orgName ] ) ) ) )

# Reading the gainLoss ancestral reconstructions.
anc_recon_table = pd.read_table( 
                  'full_proks_gainLoss_results/gainLossMP.2.00099.AncestralReconstructSankoff.txt' )

# Now also reconstructing the most likely ancestral genotypes.
def reconAncestor( anc_recon_table, node ):
    if node.name == '[N1]':
        tO = anc_recon_table.loc[ anc_recon_table['Node'] == node.name[1:-1] ]
    else:
        tO = anc_recon_table.loc[ anc_recon_table['Node'] == node.name ]

    return tO['State'].values

# Now traversing the tree and inferring ancestral states for all unmarked nodes.
for thisNode in nodes:
    try:
        thisNode.genotype
    except:
        thisNode.add_feature( 'genotype', reconAncestor( anc_recon_table, thisNode ) )

allowed_pos = np.load('full_proks_gainLoss_results/allowed_positions_genotypes.npy')
pos_to_gene_map = {i: gene_ids[i] for i in range(len(gene_ids))}
gene_to_pos_map = {value: key for key, value in pos_to_gene_map.items()}

from give_scope import giveScope

# Generating a vector of all metabolites initially provided, i.e. seeds.
seedVec = np.zeros(len(rxnMat.T))
seeds_df = pd.read_csv('seeds_from_vitkup.csv')
media = list(seeds_df['kegg_id'].values)
media = list(set(media) & set(cpds))
seedVec[[kegg_to_id[e] for e in Currency + media]] = 1

def inferIndDegree(genotype):
    kos = [pos_to_gene_map[i] for i in np.where((genotype > 0.5) * 1)[0]]
    refDF = pd.read_csv('KOREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])
    def give_number( string ):
        return int( string[-5:] )
    refDF['id'] = refDF['id'].apply( give_number )
    refDF['rid'] = refDF['rid'].apply( give_number )
    can = refDF.loc[refDF['id'].isin(kos),:]['rid'].values
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in can]
    can = list((set(can)) & set(rxns))
    avrxns = [rxn_kegg_to_id[e] for e in can]
    scopeMets = giveScope(rxnMat[avrxns], prodMat[avrxns], seedVec, sumRxnVec[avrxns])[0]
    return list(set([id_to_kegg[e] for e in set(np.where(scopeMets)[0])]) & set(Core))


# Writing all internal ancestral states.
def inferChildren(gl_node, anc_node):
    gl_node.add_feature('ind_degree', len(inferIndDegree(gl_node.genotype)))

    for gl_child, anc_child in zip(gl_node.children, anc_node.children):
        try:
            gl_child.ind_degree
        except:
            inferChildren(gl_child, anc_child)

inferChildren(root, ind_tree)

# Using first ancestral genotype inference method to calculate gains and losses.
def giveGainsAndLosses( parent, child ):
    gainGenes, lostGenes = set( ), set( )
    for indx, geneID in enumerate( gene_ids ):
        if indx in allowed_pos:
            parentProb, childProb = parent.genotype[ indx ], child.genotype[ indx ]

            # Order is present, absent, gain and loss.
            prsnProb = parentProb * childProb
            absnProb = ( 1 - parentProb ) * ( 1 - childProb )
            gainProb = ( 1 - parentProb ) * childProb
            lossProb = parentProb * ( 1 - childProb )

            # Finding most likely event by indexing.
            probList = [ prsnProb, absnProb, gainProb, lossProb ]
            mostLikelyEvent = probList.index( max( probList ) )

            # Checking if gain or loss.
            if mostLikelyEvent == 2:
                gainGenes.add( geneID )
            elif mostLikelyEvent == 3:
                lostGenes.add( geneID )
    
    return gainGenes, lostGenes

# Using first ancestral genotype inference method to calculate gains and losses.
def isGeneRetained( parent, child, geneID ):
    indx = gene_ids.index( geneID )
    parentProb, childProb = parent.genotype[ indx ], child.genotype[ indx ]

    # Order is present, absent, gain and loss.
    prsnProb = parentProb * childProb
    absnProb = ( 1 - parentProb ) * ( 1 - childProb )
    gainProb = ( 1 - parentProb ) * childProb
    lossProb = parentProb * ( 1 - childProb )

    # Finding most likely event by indexing.
    probList = [ prsnProb, absnProb, gainProb, lossProb ]
    mostLikelyEvent = probList.index( max( probList ) )

    # Checking if gain or loss.
    if mostLikelyEvent == 0:
        return True
    return False

# Getting gain to loss ratios for each group.
def giveGLRatio( transitionSet ):
    glRatios = []
    for gainGenes, lostGenes in transitionSet:
        try:
            glRatios.append( len( gainGenes ) / len( lostGenes ) )
        except:
            pass

    return glRatios

# Now traversing the tree and inferring each branch's transitions.
# delta_ind_degree = []
# nums_gains, nums_losses = [], []
# list_gains, list_losses = [], []
# gl_ratios, gl_deltas = [], []
# allGains, allLosses = [], []
# full_names = [tip_to_name_segata_dict[kegg_to_tip_name_dict[o]] for o in orgNames]
# for thisNode in nodes:
#     for thisChild in thisNode.children:
#         gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
#         allGains.append(gainGenes)
#         allLosses.append(lostGenes)
#         if (len(lostGenes) / np.count_nonzero(thisNode.genotype) >= 0.00 
#             and len(gainGenes) / np.count_nonzero(thisNode.genotype) >= 0.00):
#             delta_ind_degree.append(thisChild.ind_degree - thisNode.ind_degree)
#             nums_gains.append( len(gainGenes) )
#             nums_losses.append( len(lostGenes) )
#             list_gains.append(list(gainGenes))
#             list_losses.append(list(lostGenes))
#             gl_deltas.append(thisChild.ind_degree - thisNode.ind_degree)
#             # gl_ratios.append(len(gainGenes) / len(lostGenes))

#------------------------------------------------------------------------------

def deltaCoreAcess(sources, anc, gainedSet):
    ancCores = giveCoresProduced(sources, anc.genotype)
    addAnc = np.copy(anc.genotype)
    for p in gainedSet:
        addAnc[gene_to_pos_map[p]] = 1
    desCores = giveCoresProduced(sources, addAnc)

    return desCores, desCores.difference(ancCores)

def giveSeedVec(sources):
    # Generating a vector of all metabolites initially provided, i.e. seeds.
    seedVec = np.zeros(len(rxnMat.T))
    # media = list(sources['kegg_id'].values)
    # media = list(sources)
    media = list(sources & set(cpds))
    seedVec[[kegg_to_id[e] for e in Currency + media]] = 1
    return seedVec

def giveCoresProduced(sources, thisGenotype):
    kos = [pos_to_gene_map[i] for i in np.where((thisGenotype > 0.5) * 1)[0]]
    refDF = pd.read_csv('KOREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])
    def give_number( string ):
        return int( string[-5:] )
    refDF['id'] = refDF['id'].apply( give_number )
    refDF['rid'] = refDF['rid'].apply( give_number )
    can = refDF.loc[refDF['id'].isin(kos),:]['rid'].values
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in can]
    can = list((set(can)) & set(rxns))
    avrxns = [rxn_kegg_to_id[e] for e in can]
    scopeMets, avScopeRxn = giveScope(rxnMat[avrxns], prodMat[avrxns], giveSeedVec(sources), 
                                     sumRxnVec[avrxns])
    scopeRxns = np.nonzero(avrxns * avScopeRxn)[0]
    return set([id_to_kegg[e] for e in set(np.where(scopeMets)[0])]) & set(Core)


allByps = pd.read_csv('secreted_mets_segre.csv')
allNuts = pd.read_csv('seeds_from_vitkup.csv')
nuts = set(allNuts.iloc[:, 1].values)
byps = set(allByps.iloc[:, 1].values)
unuts = nuts.difference(byps)
ubyps = byps.difference(nuts)

# Getting core accessibility from byproducts and nutrients.
# bypUseList, nutUseList = [], []
# for thisNode in tqdm(nodes):
#     for thisChild in thisNode.children:
#         gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
#         if gainGenes:
#             bypUseList.append(deltaCoreAcess(allByps, thisNode, gainGenes))
#             nutUseList.append(deltaCoreAcess(allNuts, thisNode, gainGenes))


# Getting core accessibility from byproducts and nutrients.
bypUseList, nutUseList = [], []
for thisNode in tqdm(nodes):
    for thisChild in thisNode.children:
        gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
        if gainGenes:
            bypUseList.append(deltaCoreAcess(ubyps, thisNode, gainGenes))
            nutUseList.append(deltaCoreAcess(unuts, thisNode, gainGenes))
