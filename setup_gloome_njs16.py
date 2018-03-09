import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify

njs16_rxnDict = pickle.load( open( 'dict_njs16_rxn.dat', 'rb' ) )
njs16_ind_dict = pickle.load( open( 'dict_njs16_ind.dat', 'rb' ) )
reaction_ids = sorted( uniqify( unlistify( list( njs16_rxnDict.values() ) ) ) )
njs16_pa_dict = {}
for thisOrg in njs16_ind_dict.keys():
    njs16_pa_dict[ thisOrg ] = ''
    for thisRxnID in reaction_ids:
        if thisRxnID in njs16_rxnDict[ thisOrg ]:
            njs16_pa_dict[ thisOrg ] += '1'
        else:
            njs16_pa_dict[ thisOrg ] += '0'

njs16_names = [ '_'.join( e.split() ) for e in list( njs16_ind_dict.keys() ) ]
img_to_name_dict = pickle.load( open( 'img_to_name_dict.dat', 'rb' ) )
inv_img_to_name_dict = { value: key for key, value in img_to_name_dict.items() }

# Writing MSA presence/absence file.
with open( 'njs16_msa_file.txt', 'w' ) as msaFile:
    for tn in njs16_pa_dict:
        msaFile.write( '>' + inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() ) + '\n' + njs16_pa_dict[ tn ] + '\n' )

running_string = ''
for tn in njs16_pa_dict:
    running_string += ' \"' + inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() ) + '\",'

# Getting the pairwise matrices for phylogenetic distances, gene losses and gains.
pdMat = pd.read_csv( 'njs16_gainLoss_results/njs16_anctree_pairwise.csv' )
gsMat = pd.read_csv( 'njs16_gainLoss_results/njs16_gains_pairwise.csv' )
lsMat = pd.read_csv( 'njs16_gainLoss_results/njs16_loss_pairwise.csv' )
org_names = pdMat.loc[:, pdMat.columns[0]].values   # Here, the names are in the right order.

isIndDict = {}
for tn in njs16_ind_dict:
    orgName = inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() )
    isIndDict[ orgName ] = njs16_ind_dict[ tn ]
