import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify

NUM_ORGS = 1031


# Getting the list of important 'terminal' metabolites.
mod_to_comp_dict = pickle.load(open('module_containers/mod_to_comp_dict.dat', 'rb'))
terMets = list(set([ mod_to_comp_dict[thisMod][-1] for thisMod in mod_to_comp_dict 
                     if mod_to_comp_dict[thisMod] ]))
terMets.remove('CHPF2;')  # Spurious metabolite being removed. What the fuck is this even?


# Note that glycans have been ignored in this analysis. 
# We get a set of 228 terminal metabolites thereby.


# Getting the set of bacterial abbreviations for each organism in KEGG.
orgNames = []
with open('prok_abbr_kegg.txt', 'r') as f:
    for thisLine in f.readlines():
        orgNames.append( thisLine.strip() )


geneDict = {thisOrg : 
            set(np.genfromtxt('organism_kogenes/' + thisOrg + '.txt').astype(int))
            for thisOrg in orgNames}
gene_ids = sorted( uniqify( unlistify( list( geneDict.values() ) ) ) )
all_pa_dict = {}
for thisOrg in tqdm( orgNames ):
    all_pa_dict[ thisOrg ] = ''
    for thisGeneID in gene_ids:
        if thisGeneID in geneDict[ thisOrg ]:
            all_pa_dict[ thisOrg ] += '1'
        else:
            all_pa_dict[ thisOrg ] += '0'

# Filtering out bad indices, that have only 1 value for all genomes.
bad_indices = []
for i in range(len(gene_ids)):
    if len( set( [ int(all_pa_dict[tn][i]) for tn in all_pa_dict ] ) ) == 1:
        bad_indices.append( i )

# Now getting the good indices and modifying the pa dict.
good_indices = sorted([ x for x in list( range( len( gene_ids ) ) ) if x not in bad_indices ])

all_pa_dict = {}
for thisOrg in tqdm( orgNames ):
    all_pa_dict[ thisOrg ] = ''
    for thisGeneID in good_indices:
        if thisGeneID in geneDict[ thisOrg ]:
            all_pa_dict[ thisOrg ] += '1'
        else:
            all_pa_dict[ thisOrg ] += '0'

# Loading all the mapping files generated hopelessly.
kegg_to_tip_name_dict = pickle.load(open('kegg_to_tip_name_dict.dat', 'rb'))
tip_to_name_segata_dict = pickle.load(open('tip_to_name_segata_dict.dat', 'rb'))
# tip_to_name_segata_dict[kegg_to_tip_name_dict[thisOrg]]

# Writing MSA presence/absence file.
# num_failed = 0
# with open( 'all_msa_file.txt', 'w' ) as msaFile:
#     for thisOrg in orgNames:
#         try:
#             msaFile.write( '>' + tip_to_name_segata_dict[kegg_to_tip_name_dict[thisOrg]] + '\n' + all_pa_dict[ thisOrg ] + '\n' )
#         except:
#             num_failed += 1
#             continue

# running_string = ''
# for thisOrg in orgNames:
#     try:
#         running_string += ' \"' + tip_to_name_segata_dict[kegg_to_tip_name_dict[thisOrg]] + '\",'
#     except:
#         continue

# isIndDict = {}
# num_failed = 0
# for thisOrg in FULLPROTRAITSisThereDict:
#     try:
#         orgName = inv_img_to_name_dict[ '_'.join( orgNameDict[thisOrg].split() ) ] + '_' + '_'.join( orgNameDict[thisOrg].split() )
#         isIndDict[ orgName ] = FULLPROTRAITSisThereDict[ thisOrg ]
#     except:
#         num_failed += 1

# genotype_dict = {}
# num_failed = 0
# for thisOrg in all_pa_dict:
#     try:
#         orgName = inv_img_to_name_dict[ '_'.join( thisOrg.split() ) ] + '_' + '_'.join( thisOrg.split() )
#         genotype_dict[ orgName ] = all_pa_dict[ thisOrg ]
#     except:
#         num_failed += 1
