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
delta_ind_degree = []
nums_gains, nums_losses = [], []
list_gains, list_losses = [], []
gl_ratios, gl_deltas = [], []
allGains, allLosses = [], []
full_names = [tip_to_name_segata_dict[kegg_to_tip_name_dict[o]] for o in orgNames]
for thisNode in nodes:
    for thisChild in thisNode.children:
        gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
        allGains.append(gainGenes)
        allLosses.append(lostGenes)
        if (len(lostGenes) / np.count_nonzero(thisNode.genotype) >= 0.00 
            and len(gainGenes) / np.count_nonzero(thisNode.genotype) >= 0.00):
            delta_ind_degree.append(thisChild.ind_degree - thisNode.ind_degree)
            nums_gains.append( len(gainGenes) )
            nums_losses.append( len(lostGenes) )
            list_gains.append(list(gainGenes))
            list_losses.append(list(lostGenes))
            gl_deltas.append(thisChild.ind_degree - thisNode.ind_degree)
            # gl_ratios.append(len(gainGenes) / len(lostGenes))
            





import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import statsmodels.api as sm
fig, ax = plt.subplots(1)
ax.scatter(delta_ind_degree, nums_gains, c='g')
lowess = sm.nonparametric.lowess(nums_gains, delta_ind_degree, frac=.6)
lowess_x = list(zip(*lowess))[0]
lowess_y = list(zip(*lowess))[1]
f = interp1d(lowess_x, lowess_y, bounds_error=False)
xnew = np.linspace(min(gl_deltas), max(gl_deltas), 1000)
ynew = f(xnew)
ax.plot(xnew, ynew, c='grey', lw=3)
ax.set_ylim(0.0, max(nums_losses)+10)
plt.savefig('dists/mapping_anc_inferred_gains_delta_i.svg')
plt.show()

fig, ax = plt.subplots(1)
ax.scatter(delta_ind_degree, nums_losses, c='r')
lowess = sm.nonparametric.lowess(nums_losses, delta_ind_degree, frac=.6)
lowess_x = list(zip(*lowess))[0]
lowess_y = list(zip(*lowess))[1]
f = interp1d(lowess_x, lowess_y, bounds_error=False)
xnew = np.linspace(min(gl_deltas), max(gl_deltas), 1000)
ynew = f(xnew)
ax.plot(xnew, ynew, c='grey', lw=3)
ax.set_ylim(0.0, max(nums_losses)+10)
plt.savefig('dists/mapping_anc_inferred_losses_delta_i.svg')
plt.show()

fig, ax = plt.subplots(1)
ax.scatter(gl_deltas, gl_ratios, c='r')
lowess = sm.nonparametric.lowess(gl_ratios, gl_deltas, frac=.7)
lowess_x = list(zip(*lowess))[0]
lowess_y = list(zip(*lowess))[1]
f = interp1d(lowess_x, lowess_y, bounds_error=False)
xnew = np.linspace(min(gl_deltas), max(gl_deltas), 1000)
ynew = f(xnew)
ax.plot(xnew, ynew, c='grey', lw=3)
plt.savefig('dists/mapping_anc_inferred_glRatios_full.svg')
plt.show()


# # Saved ind_to_dep_list as a pickle object.
# # Going down subsequent branches of transiting origins and measure gain/loss retention.
# num_retain_trans_list, num_retain_control_list = [], []
# num_gains_trans_list, num_gains_control_list = [], []
# prob_retain_trans, prob_retain_control = [], []
# for tI, thisTrans in enumerate( ind_to_dep_list ):
#     # Getting the gained genes list for this transition.
#     gainsSet = ind_to_dep_GLS[ tI ][ 0 ]
#     if not gainsSet:
#         continue

#     if thisTrans[1].children:
#         thisChild = (thisTrans[1].children)[0]
#         retainedSet = set( [ gID 
#                              for gID in gainsSet
#                              if isGeneRetained( thisTrans[1], thisChild, gID ) ] )
#         num_retain_trans = len( retainedSet )
#         num_retain_trans_list.append( num_retain_trans )
#         num_gains_trans_list.append( len( gainsSet ) )
#         prob_retain_trans.append( num_retain_trans / len( gainsSet ) )

# for tI, thisTrans in enumerate( ind_to_ind_list ):
#     # Getting the gained genes list for this transition.
#     gainsSet = ind_to_ind_GLS[ tI ][ 0 ]
#     gainsInDep = bool( ind_to_dep_GLS[ tI ][ 0 ] )
#     if not gainsSet or not gainsInDep:
#         continue

#     if thisTrans[1].children:
#         thisChild = (thisTrans[1].children)[0]
#         retainedSet = set( [ gID 
#                              for gID in gainsSet
#                              if isGeneRetained( thisTrans[1], thisChild, gID ) ] )
#         num_retain_control = len( retainedSet )
#         num_retain_control_list.append( num_retain_control )
#         num_gains_control_list.append( len( gainsSet ) )
#         prob_retain_control.append( num_retain_control / len( gainsSet ) )

# # Getting the pathway representation for the gains.
# isModuleComplete = []
# for theseGains, theseLosses in ind_to_ind_GLS:
#     pr_expr_dict = pickle.load( open( 'module_exprs/pr_expr_dict.dat', 'rb' ) )
#     isModuleComplete.append( np.array( [] ) )
#     for thisMod in pr_expr_dict:
#         for x in re.findall( r'K\d{5}', pr_expr_dict[ thisMod ] ):
#             pr_expr_dict[ thisMod ] = pr_expr_dict[ thisMod ].replace( 
#                                       x, str( int( x[ 1: ] ) in theseGains ) )
#         isModuleComplete[ -1 ] = np.append( isModuleComplete[ -1 ], 
#                                             eval( pr_expr_dict[ thisMod ] ) )

# control_modules_complete = unlistify( list(np.where(isModuleComplete[x])[0]) for x in range(len(isModuleComplete)) )

# # Getting the pathway representation for the gains.
# isModuleComplete = []
# for theseGains, theseLosses in ind_to_dep_GLS:
#     pr_expr_dict = pickle.load( open( 'module_exprs/pr_expr_dict.dat', 'rb' ) )
#     isModuleComplete.append( np.array( [] ) )
#     for thisMod in pr_expr_dict:
#         for x in re.findall( r'K\d{5}', pr_expr_dict[ thisMod ] ):
#             pr_expr_dict[ thisMod ] = pr_expr_dict[ thisMod ].replace( 
#                                       x, str( int( x[ 1: ] ) in theseGains ) )
#         isModuleComplete[ -1 ] = np.append( isModuleComplete[ -1 ], 
#                                             eval( pr_expr_dict[ thisMod ] ) )

# test_modules_complete = unlistify( list(np.where(isModuleComplete[x])[0]) 
#                                       for x in range(len(isModuleComplete)) )

# # Getting the number of gains and losses in test and control branches.
# control_num_gains = [ len( ind_to_ind_GLS[x][0] ) for x in range(len(ind_to_ind_GLS)) ]
# test_num_gains = [ len( ind_to_dep_GLS[x][0] ) for x in range(len(ind_to_dep_GLS)) ]
# control_num_losses = [ len( ind_to_ind_GLS[x][1] ) for x in range(len(ind_to_ind_GLS)) ]
# test_num_losses = [ len( ind_to_dep_GLS[x][1] ) for x in range(len(ind_to_dep_GLS)) ]

def inferByps(genotype):
    allByps = pd.read_csv('secreted_mets_segre.csv')
    byps = set(allByps.iloc[:, 1].values)
    ubyps = byps.difference(nuts)

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
    return list(set([id_to_kegg[e] for e in set(np.where(scopeMets)[0])]))

max_diversity = 10
bypsAtDiv = { e : set([]) for e in range(2, max_diversity + 1)}
for currDiv in range(2, max_diversity + 1):
    chosenNodes = random.sample(nodes, currDiv)
    for thisNode in chosenNodes:
        thisNodeByps = set(inferByps( thisNode.genotype ))
        bypsAtDiv[currDiv] = bypsAtDiv[currDiv].union(thisNodeByps)

sns.lineplot(np.array(list(range(2,max_diversity + 1))), list(map(len, list(bypsAtDiv.values()))))
plt.show()