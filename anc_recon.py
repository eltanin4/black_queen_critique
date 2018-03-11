from ete3 import Tree
import numpy as np
from setup_gloome_njs16 import isIndDict, genotype_dict
import pandas as pd

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
for orgName in isIndDict:
        (ancTree_njs16&orgName).add_feature( 'genotype', genotype_dict[ orgName ] )

# Reading the gainLoss ancestral reconstructions.
anc_recon_table = pd.read_table( 
                  'njs16_gainLoss_results/RESULTS/AncestralReconstructPosterior.txt' )

# Now also reconstructing the most likely ancestral genotypes.
def reconAncestor( anc_recon_table, node ):
    if node.name == '[N1]':
        tO = anc_recon_table.loc[ anc_recon_table['Node'] == node.name[1:-1] ]
    else:
        tO = anc_recon_table.loc[ anc_recon_table['Node'] == node.name ]

    return ( tO[ 'Prob' ].values > 0.5 ) * 1
