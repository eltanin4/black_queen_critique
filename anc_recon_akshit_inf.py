from ete3 import Tree
import numpy as np
from setup_gloome_param_files import isIndDict, genotype_dict, good_indices
import pandas as pd
from uniqify import uniqify
from unlistify import unlistify
import pickle

# Getting the tree with internal nodes from gainLoss' ancestral reconstruction.
# Internal nodes are labeled with 'Nx' where x is a number, the root being
# '[N1]' and then the numbers increase.
ancTree_njs16 = Tree( 'full_proks_gainLoss_results/TheTree.INodes.ph', format=1 )

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
# for orgName in isIndDict:
#         (ancTree_njs16&orgName).add_feature( 'genotype', genotype_dict[ orgName ] )

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

njs16_rxnDict = pickle.load( open( 'dict_njs16_rxn.dat', 'rb' ) )
gene_ids = sorted( uniqify( unlistify( list( njs16_rxnDict.values() ) ) ) )
gene_ids = list( np.array( gene_ids )[ good_indices ] )

# Using first ancestral genotype inference method to calculate gains and losses.
def giveGainsAndLosses( parent, child ):
    # Converting the binary strings to arrays.
    # parentState = np.array( list( map( int, parent.genotype ) ) )
    # childState = np.array( list( map( int, child.genotype ) ) )

    # Getting the array of reaction IDs present in both parent and child.
    parentRxns = np.nonzero( np.multiply( parent.genotype, gene_ids ) )[0]
    childRxns = np.nonzero( np.multiply( child.genotype, gene_ids ) )[0]
    
    commonRxns = set( parentRxns ) & set( childRxns )
    lostRxns = set( parentRxns ).difference( commonRxns )
    gainRxns = set( childRxns ).difference( commonRxns )

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

# Now to try the alternate 'probability' based way to infer gains and losses 
# (as in Press et. al., Gen. Res. 2016).
# Getting gain to loss ratios for each group.
def giveGLRatio( transitionSet ):
    glRatios = []
    for gainRxns, lostRxns in transitionSet:
        if len( lostRxns ) and len( gainRxns ):
            glRatios.append( len( gainRxns ) / len( lostRxns ) )
        else:
            continue

    return glRatios

# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(1)
# plt.hist( giveGLRatio( ind_to_dep_list ), color='dodgerblue' )
# plt.hist( giveGLRatio( dep_to_dep_list ), color='mediumseagreen' )
# plt.hist( giveGLRatio( dep_to_dep_list ), color='firebrick' )
# plt.show()
