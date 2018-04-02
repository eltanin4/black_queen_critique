from ete3 import Tree
import numpy as np
from setup_gloome_param_files import isIndDict, genotype_dict, geneDict, good_indices
import pandas as pd
from uniqify import uniqify
from unlistify import unlistify
from copy import copy
import pickle
import re

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
# NOTE: genotype dicts are now probability vectors.
# for orgName in isIndDict:
#         (ancTree_njs16&orgName).add_feature( 'genotype', 
#                                 np.array( list( map( int, genotype_dict[ orgName ] ) ) ) )

# Reading the gainLoss ancestral reconstructions.
anc_recon_table = pd.read_table( 
                  'full_proks_gainLoss_results/AncestralReconstructPosterior.txt' )

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


gene_ids = sorted( uniqify( unlistify( list( geneDict.values() ) ) ) )
gene_ids = list( np.array( gene_ids )[ good_indices ] )

# Using first ancestral genotype inference method to calculate gains and losses.
def giveGainsAndLosses( parent, child ):
    gainGenes, lostGenes = set( ), set( )
    for indx, geneID in enumerate( gene_ids ):
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
ind_to_dep_GLS, ind_to_ind_GLS, dep_to_dep_GLS = [], [], []
ind_to_dep_list, ind_to_ind_list, dep_to_dep_list = [], [], []
for thisNode in nodes:
    for thisChild in thisNode.children:
        if thisNode.isInd and not thisChild.isInd:
            print( thisNode.name, thisChild.name )
            ind_to_dep_list.append( [ thisNode, thisChild ] )
            ind_to_dep_GLS.append( giveGainsAndLosses( thisNode, thisChild ) )
            
            # Test if control branch exists.
            for otherChild in thisNode.children:
                if otherChild != thisChild and otherChild.isInd:
                    ind_to_ind_list.append( [ thisNode, otherChild ] )
                    ind_to_ind_GLS.append( giveGainsAndLosses( thisNode, otherChild ) )

# Saved ind_to_dep_list as a pickle object.
# Going down subsequent branches of transiting origins and measure gain/loss retention.
num_retain_trans_list, num_retain_control_list = [], []
num_gains_trans_list, num_gains_control_list = [], []
prob_retain_trans, prob_retain_control = [], []
for tI, thisTrans in enumerate( ind_to_dep_list ):
    # Getting the gained genes list for this transition.
    gainsSet = ind_to_dep_GLS[ tI ][ 0 ]
    if not gainsSet:
        continue

    if thisTrans[1].children:
        thisChild = (thisTrans[1].children)[0]
        retainedSet = set( [ gID 
                             for gID in gainsSet
                             if isGeneRetained( thisTrans[1], thisChild, gID ) ] )
        num_retain_trans = len( retainedSet )
        num_retain_trans_list.append( num_retain_trans )
        num_gains_trans_list.append( len( gainsSet ) )
        prob_retain_trans.append( num_retain_trans / len( gainsSet ) )

for tI, thisTrans in enumerate( ind_to_ind_list ):
    # Getting the gained genes list for this transition.
    gainsSet = ind_to_ind_GLS[ tI ][ 0 ]
    gainsInDep = bool( ind_to_dep_GLS[ tI ][ 0 ] )
    if not gainsSet or not gainsInDep:
        continue

    if thisTrans[1].children:
        thisChild = (thisTrans[1].children)[0]
        retainedSet = set( [ gID 
                             for gID in gainsSet
                             if isGeneRetained( thisTrans[1], thisChild, gID ) ] )
        num_retain_control = len( retainedSet )
        num_retain_control_list.append( num_retain_control )
        num_gains_control_list.append( len( gainsSet ) )
        prob_retain_control.append( num_retain_control / len( gainsSet ) )

# Getting the pathway representation for the gains.
isModuleComplete = []
for theseGains, theseLosses in ind_to_ind_GLS:
    pr_expr_dict = pickle.load( open( 'module_exprs/pr_expr_dict.dat', 'rb' ) )
    isModuleComplete.append( np.array( [] ) )
    for thisMod in pr_expr_dict:
        for x in re.findall( r'K\d{5}', pr_expr_dict[ thisMod ] ):
            pr_expr_dict[ thisMod ] = pr_expr_dict[ thisMod ].replace( 
                                      x, str( int( x[ 1: ] ) in theseGains ) )
        isModuleComplete[ -1 ] = np.append( isModuleComplete[ -1 ], 
                                            eval( pr_expr_dict[ thisMod ] ) )

control_modules_complete = unlistify( list(np.where(isModuleComplete[x])[0]) for x in range(len(isModuleComplete)) )

# Getting the pathway representation for the gains.
isModuleComplete = []
for theseGains, theseLosses in ind_to_dep_GLS:
    pr_expr_dict = pickle.load( open( 'module_exprs/pr_expr_dict.dat', 'rb' ) )
    isModuleComplete.append( np.array( [] ) )
    for thisMod in pr_expr_dict:
        for x in re.findall( r'K\d{5}', pr_expr_dict[ thisMod ] ):
            pr_expr_dict[ thisMod ] = pr_expr_dict[ thisMod ].replace( 
                                      x, str( int( x[ 1: ] ) in theseGains ) )
        isModuleComplete[ -1 ] = np.append( isModuleComplete[ -1 ], 
                                            eval( pr_expr_dict[ thisMod ] ) )

test_modules_complete = unlistify( list(np.where(isModuleComplete[x])[0]) 
                                      for x in range(len(isModuleComplete)) )

# Getting the number of gains and losses in test and control branches.
control_num_gains = np.mean( [ len( ind_to_ind_GLS[x][0] ) for x in range(len(ind_to_ind_GLS)) ] )
test_num_gains = np.mean( [ len( ind_to_dep_GLS[x][0] ) for x in range(len(ind_to_dep_GLS)) ] )
control_num_losses = np.mean( [ len( ind_to_ind_GLS[x][1] ) for x in range(len(ind_to_ind_GLS)) ] )
test_num_losses = np.mean( [ len( ind_to_dep_GLS[x][1] ) for x in range(len(ind_to_dep_GLS)) ] )
