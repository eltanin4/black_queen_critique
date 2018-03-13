import pickle

gene_ids = sorted( uniqify( unlistify( list( geneDict.values() ) ) ) )
all_pa_dict = {}
for thisOrg in tqdm( FULLPROTRAITSisThereDict.keys() ):
    all_pa_dict[ orgNameDict[ thisOrg ] ] = ''
    for thisGeneID in gene_ids:
        if thisGeneID in geneDict[ thisOrg ]:
            all_pa_dict[ orgNameDict[ thisOrg ] ] += '1'
        else:
            all_pa_dict[ orgNameDict[ thisOrg ] ] += '0'

# Filtering out bad indices, that have only 1 value for all genomes.
bad_indices = []
for i in range(len(gene_ids)):
    if len( set( [ int(all_pa_dict[tn][i]) for tn in all_pa_dict ] ) ) == 1:
        bad_indices.append( i )

# Now getting the good indices and modifying the pa dict.
good_indices = [ x for x in list( range( len( gene_ids ) ) ) if x not in bad_indices ]

for tn in all_pa_dict:
    all_pa_dict[tn] = ''.join(np.array(list(map(int, 
                           all_pa_dict[tn])))[good_indices].astype(str))

# all_names = [ '_'.join( e.split() ) for e in list( all_ind_dict.keys() ) ]
inv_img_to_name_dict = pickle.load( open( 'inv_img_to_name_dict.dat', 'rb' ) )

# Writing MSA presence/absence file.
# num_failed = 0
# with open( 'all_msa_file.txt', 'w' ) as msaFile:
#     for tn in all_pa_dict:
#         try:
#             msaFile.write( '>' + inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() ) + '\n' + all_pa_dict[ tn ] + '\n' )
#         except:
#             num_failed += 1
#             continue

running_string = ''
for tn in all_pa_dict:
    try:
        running_string += ' \"' + inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() ) + '\",'
    except:
        continue
