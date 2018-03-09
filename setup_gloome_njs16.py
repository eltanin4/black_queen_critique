import pickle

njs16_rxnDict = pickle.load( open( 'dict_njs16_rxn.dat', 'rb' ) )
njs16_ind_Dict = pickle.load( open( 'dict_njs16_ind.dat', 'rb' ) )
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
