import pickle

# 4,123 reactions
reaction_ids = sorted( uniqify( unlistify( list( rxnDict.values() ) ) ) )
all_pa_dict = {}
for thisOrg in tqdm( FULLPROTRAITSisThereDict.keys() ):
    all_pa_dict[ orgNameDict[ thisOrg ] ] = ''
    for thisRxnID in reaction_ids:
        if thisRxnID in rxnDict[ thisOrg ]:
            all_pa_dict[ orgNameDict[ thisOrg ] ] += '1'
        else:
            all_pa_dict[ orgNameDict[ thisOrg ] ] += '0'

# njs16_names = [ '_'.join( e.split() ) for e in list( njs16_ind_dict.keys() ) ]
img_to_name_dict = pickle.load( open( 'img_to_name_dict.dat', 'rb' ) )
inv_img_to_name_dict = { value: key for key, value in img_to_name_dict.items() }

# Writing MSA presence/absence file.
num_failed = 0
with open( 'all_msa_file.txt', 'w' ) as msaFile:
    for tn in all_pa_dict:
        try:
            msaFile.write( '>' + inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() ) + '\n' + all_pa_dict[ tn ] + '\n' )
        except:
            num_failed += 1
            continue

running_string = ''
for tn in all_pa_dict:
    try:
        running_string += ' \"' + inv_img_to_name_dict[ '_'.join( tn.split() ) ] + '_' + '_'.join( tn.split() ) + '\",'
    except:
        continue
