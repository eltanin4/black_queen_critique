from ete3 import Tree

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

markNode( root )
