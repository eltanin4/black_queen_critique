from ete3 import Tree
import numpy as np
from setup_gloome_njs16 import isIndDict, genotype_dict
import pandas as pd
from uniqify import uniqify
from unlistify import unlistify
from copy import copy
import pickle

# Getting the tree with internal nodes from gainLoss' ancestral reconstruction.
# Internal nodes are labeled with 'Nx' where x is a number, the root being
# '[N1]' and then the numbers increase.
ancTree_njs16 = Tree( 'njs16_gainLoss_results/RESULTS/TheTree.INodes.ph', format=1 )

# To print the tree.
print(ancTree_njs16.get_ascii(show_internal=True))

# Getting the root of the tree.
root = ancTree_njs16&'[N1]'
nodes = list( root.traverse() )

# Getting the original traits in a trait dict.
for orgName in isIndDict:
        (ancTree_njs16&orgName).add_feature( 'isInd', isIndDict[ orgName ] )

# Writing the markNode recursive function.
def markNode( tree, node ):
    children = node.children

    # Checking if all children marked.
    for thisChild in children:
        try:
            thisChild.isInd
        except:
            markNode( tree, thisChild )

    # When you get here, all your children should be marked.
    if np.array( [ not thisChild.isInd for thisChild in children ] ).all():
        node.add_feature( 'isInd', 0 )
    else:
        node.add_feature( 'isInd', 1 )

markNode( ancTree_njs16, root )

# Getting the original genotypes in a genotype dict.
# NOTE: genotype dicts are now probability vectors.
# for orgName in isIndDict:
#         (ancTree_njs16&orgName).add_feature( 'genotype', 
#                                 np.array( list( map( int, genotype_dict[ orgName ] ) ) ) )

# Reading the gainLoss ancestral reconstructions.
anc_recon_table = pd.read_table( 
                  'njs16_gainLoss_results/RESULTS/AncestralReconstructPosterior.txt' )

# Now to try the alternate 'probability' based way to infer gains and losses 
# (as in Press et. al., Gen. Res. 2016).
def reconAncestor( anc_recon_table, node ):
    if node.name == '[N1]':
        tO = anc_recon_table.loc[ anc_recon_table['Node'] == node.name[1:-1] ]
    else:
        tO = anc_recon_table.loc[ anc_recon_table['Node'] == node.name ]

    return tO['Prob'].values

# Now traversing the tree and inferring ancestral states for all unmarked nodes.
for thisNode in nodes:
    try:
        thisNode.genotype
    except:
        thisNode.add_feature( 'genotype', reconAncestor( anc_recon_table, thisNode ) )


njs16_rxnDict = pickle.load( open( 'dict_njs16_rxn.dat', 'rb' ) )
reaction_ids = sorted( uniqify( unlistify( list( njs16_rxnDict.values() ) ) ) )

# Using first ancestral genotype inference method to calculate gains and losses.
def giveGainsAndLosses( parent, child ):
    gainRxns, lostRxns = set( ), set( )
    for indx, rxnID in enumerate( reaction_ids ):
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
            gainRxns.add( rxnID )
        elif mostLikelyEvent == 3:
            lostRxns.add( rxnID )
    
    return gainRxns, lostRxns

# Now traversing the tree and inferring each branch's transitions.
ind_to_dep_list, ind_to_ind_list, dep_to_dep_list = [], [], []
for thisNode in nodes:
    for thisChild in thisNode.children:
        if thisNode.isInd and thisChild.isInd:
            ind_to_ind_list.append( giveGainsAndLosses( thisNode, thisChild ) )
        elif not thisNode.isInd and not thisChild.isInd:
            dep_to_dep_list.append( giveGainsAndLosses( thisNode, thisChild ) )
        elif thisNode.isInd and not thisChild.isInd:
            print( thisNode.name, thisChild.name )
            ind_to_dep_list.append( giveGainsAndLosses( thisNode, thisChild ) )

# Getting gain to loss ratios for each group.
def giveGLRatio( transitionSet ):
    glRatios = []
    for gainRxns, lostRxns in transitionSet:
        try:
            glRatios.append( len( gainRxns ) / len( lostRxns ) )
        except:
            pass

    return glRatios

# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(1)
# plt.hist( giveGLRatio( ind_to_dep_list ), color='dodgerblue' )
# plt.hist( giveGLRatio( dep_to_dep_list ), color='mediumseagreen' )
# plt.hist( giveGLRatio( dep_to_dep_list ), color='firebrick' )
# plt.show()
