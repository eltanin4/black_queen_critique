import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify
NUM_ORGS = 301

# Getting the names of the abbrevations from KEGG for KOMODO comparison.
orgNameDict = {}
with open('prok_kegg_abbr_with_names_dups_removed.txt', 'r') as f:
    for thisOrg in range(NUM_ORGS):              
        thisLine = f.readline()    
        orgNameDict[ thisLine[:3] ] = thisLine[4:-1]

prokNames = list( orgNameDict.keys() ) 

# Getting sizes for each organism.
geneDict = {}
for thisOrg in prokNames:
    geneDict[ thisOrg ] = list( np.genfromtxt( 'organism_kogenes/' + thisOrg + '.txt' ) )

# # March 2, 2018
# # Getting KOMODO growers again.
# import pandas as pd
org_to_media_df = pd.read_html('komodo_org_to_media.html')[1]

# March 5, 2018
for thisInd, thisOrg in enumerate( org_to_media_df[2][1:] ):
    try:
        org_to_media_df[2][thisInd] = ' '.join( thisOrg.split( )[ : 2 ] )
    except:
        pass

# Computing whether a medium is listed for a species or not.
KOMODOisThereDict = {}
for thisOrg in prokNames:
    KOMODOisThereDict[ thisOrg ] = bool( org_to_media_df[2].str.contains( 
                                   orgNameDict[ thisOrg ] ).any() )

# Getting the normalized and manually checked ProTraits growers.
invOrgNameDict = { value: key for key, value in orgNameDict.items() }
FULLPROTRAITSisThereDict = {}
proCultDF = pd.read_csv('filtered_ProTraits_binaryIntegratedPr0.95.csv')
for idx, thisRow in proCultDF.iterrows():
    try:
        thisOrgName = invOrgNameDict[ ' '.join( thisRow[1].split( )[ : 2 ] ) ]
    except:
        continue

    if bool( thisRow[-1] == True ):
        FULLPROTRAITSisThereDict[ thisOrgName ] = 1
    else:
        try:
            FULLPROTRAITSisThereDict[ thisOrgName ]
        except:
            FULLPROTRAITSisThereDict[ thisOrgName ] = 0

# Getting sizes for each organism.
geneDict = {}
for thisOrg in FULLPROTRAITSisThereDict.keys():
    geneDict[ thisOrg ] = list( np.genfromtxt( 'organism_kogenes/' + thisOrg + '.txt' ) )

# Comparing sizes.
sizes_ind = [ len( geneDict[ thisOrgName ] ) 
              for thisOrgName in geneDict.keys() 
              if FULLPROTRAITSisThereDict[ thisOrgName ] ] 
sizes_dep = [ len( geneDict[ thisOrgName ] ) 
              for thisOrgName in geneDict.keys() 
              if not FULLPROTRAITSisThereDict[ thisOrgName ] ] 
