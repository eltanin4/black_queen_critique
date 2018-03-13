import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify

njs16_geneDict = pickle.load( open( 'dict_njs16_gene.dat', 'rb' ) )
njs16_ind_dict = pickle.load( open( 'dict_njs16_ind.dat', 'rb' ) )
gene_ids = sorted( uniqify( unlistify( list( njs16_geneDict.values() ) ) ) )
njs16_pa_dict = {}
for thisOrg in njs16_ind_dict.keys():
    njs16_pa_dict[ thisOrg ] = ''
    for thisGeneID in gene_ids:
        if thisGeneID in njs16_geneDict[ thisOrg ]:
            njs16_pa_dict[ thisOrg ] += '1'
        else:
            njs16_pa_dict[ thisOrg ] += '0'

# Filtering out bad indices, that have only 1 value for all genomes.
bad_indices = []
for i in range(len(gene_ids)):
    if len( set( [ int(njs16_pa_dict[tn][i]) for tn in njs16_pa_dict ] ) ) == 1:
        bad_indices.append( i )

# Now getting the good indices and modifying the pa dict.
good_indices = [ x for x in list( range( len( gene_ids ) ) ) if x not in bad_indices ]

for tn in njs16_pa_dict:
    njs16_pa_dict[tn] = ''.join(np.array(list(map(int, 
                           njs16_pa_dict[tn])))[good_indices].astype(str))

njs16_names = [ '_'.join( e.split() ) for e in list( njs16_ind_dict.keys() ) ]
inv_img_to_name_dict = pickle.load( open( 'njs16_inv_img_to_name_dict.dat', 'rb' ) )

# Writing MSA presence/absence file.
# with open( 'njs16_msa_file.txt', 'w' ) as msaFile:
#     for tn in njs16_pa_dict:
#         msaFile.write( '>' + inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() ) + '\n' + njs16_pa_dict[ tn ] + '\n' )

running_string = ''
for tn in njs16_pa_dict:
    running_string += ' \"' + inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() ) + '\",'

# Getting the pairwise matrices for phylogenetic distances, gene losses and gains.
# pdMat = pd.read_csv( 'njs16_gainLoss_results/njs16_anctree_pairwise.csv' )
# gsMat = pd.read_csv( 'njs16_gainLoss_results/njs16_gains_pairwise.csv' )
# lsMat = pd.read_csv( 'njs16_gainLoss_results/njs16_loss_pairwise.csv' )
# org_names = pdMat.loc[:, pdMat.columns[0]].values   # Here, the names are in the right order.

isIndDict = {}
for tn in njs16_ind_dict:
    orgName = inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() )
    isIndDict[ orgName ] = njs16_ind_dict[ tn ]

genotype_dict = {}
for tn in njs16_pa_dict:
    orgName = inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() )
    genotype_dict[ orgName ] = njs16_pa_dict[ tn ]
