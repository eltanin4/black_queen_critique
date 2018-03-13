import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify
NUM_ORGS = 301

 # Getting all organism seed strings
# orgSeedList = []
#  with open('seeds.txt', 'r') as f:
#      for thisOrg in range(NUM_ORGS):              
#          thisLine = f.readline()    
#          orgSeedList.append(thisLine)

# def find_between_r(s, first, last):
#     try:
#         start = s.rindex(first) + len(first)
#         end = s.rindex(last, start)
#         return s[start:end]
#     except ValueError:
#         return ""

# # Converting all strings to list of list of numbers, each ID one of the KEGG seed IDs.
# orgSeedList = [ [ int(e) for e in find_between_r( orgSeed, ':', '\n').split(' ')[1:-1] ] for orgSeed in orgSeedList ]

# pickle.dump( orgSeedList, open( 'org_seed_list.dat', 'wb' ) )

# Getting the names of the abbrevations from KEGG for KOMODO comparison.
orgNameDict = {}
with open('prok_kegg_abbr_with_names_dups_removed.txt', 'r') as f:
    for thisOrg in range(NUM_ORGS):              
        thisLine = f.readline()    
        orgNameDict[ thisLine[:3] ] = thisLine[4:-1]

prokNames = list( orgNameDict.keys() ) 

# Getting sizes for each organism.
rxnDict = {}
for thisOrg in prokNames:
    rxnDict[ thisOrg ] = list( np.genfromtxt( 'organism_reactions/' + thisOrg + '.txt' ) )

# Get the dependency score distribution and plotting it with some standard organism values.
# depScores = list( pickle.load( open( 'dict_org_dep_scores.dat', 'rb' ) ).values() )
# normedDepWeights = np.ones_like( depScores ) / len( depScores )
# fig, ax = plt.subplots(1)
# plt.hist( depScores, bins=40, histtype='stepfilled', weights=normedDepWeights )
# ax.axvline( depScoreDict[ 'eco' ], color='black', linewidth=2.5)
# plt.savefig( 'dists/init_dep_score_dict_with_eco.svg' )
# plt.show()

# March 1, 2018
# MEM_CUTOFF = 0.95
# seedMembershipDict = {}
# all_seeds = uniqify( unlistify( orgSeedList ) )
# for thisSeed in all_seeds:
#     seedMembershipDict[ thisSeed ] = sum( [ bool( thisSeed in orgSeedDict[ thisOrg ] ) 
#                                           for thisOrg in prokNames ] )

# # Seed membership histogram is bimodal (U-shaped), so extracting top 10% of seed compounds.
# PERCENT_CHOSEN = 10
# listSeedMembership = list( seedMembershipDict.values() )
# coreSeeds = sorted( range( len( listSeedMembership ) ), 
#                     key=lambda i: listSeedMembership[ i ], reverse=True )[ : 
#                     int( PERCENT_CHOSEN / 100 * len( prokNames ) ) ]

# # Now removing all core seeds to make a non-core seed dict for all organisms.
# orgSeedDict_CORE = {}
# for thisOrg in prokNames:
#     orgSeedDict_CORE[ thisOrg ] = [ thisSeed for thisSeed in orgSeedDict[ thisOrg ] 
#                                     if thisSeed not in coreSeeds]

# Now getting the "core seed normalized" dependency scores.
# depScoreDict_CORE = {}
# for thisOrg in prokNames:
#     depScoreDict_CORE[ thisOrg ] = ( len( orgSeedDict_CORE[ thisOrg ] ) / 
#                                      len( orgRxnDict[ thisOrg ] ) )

# # Now plotting all core normalized dependency scores.
# depScores_CORE = list( depScoreDict_CORE.values() )
# depScores_COREWeights = np.ones_like( np.array( depScores_CORE ) ) / len( depScores_CORE )
# plt.hist( depScores_CORE, histtype='stepfilled', bins=40, weights=depScores_COREWeights )
# plt.show()

# Getting ProTraits grower list.
# proGrowersList = []
# with open('proTraits_growers.txt', 'r') as f:
#     for thisLine in f:
#         proGrowersList.append( thisLine[:-1] )

# Checking against pro-traits growers on many many media.
# PROTRAITSisThereDict = {}
# for thisOrg in prokNames:
#     PROTRAITSisThereDict[ thisOrg ] = bool( orgNameDict[ thisOrg ] in proGrowersList )

# Checking dependency scores for all ProTraits predicted organisms.
# depScoreDict = pickle.load( open( 'dict_org_dep_scores.dat', 'rb' ) )
# depScores = list( depScoreDict.values() )
# indDepScoreDict = { thisOrg: depScoreDict[ thisOrg ] 
#                     for thisOrg in prokNames 
#                     if  PROTRAITSisThereDict[ thisOrg ] }

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

# # Double checking if numpy is fucking up here.
# # Yup, numpy fucked up. The concordance score now looks like 61%.
# concordance_score = 0.0
# for thisOrg in prokNames:
#     if KOMODOisThereDict[ thisOrg ] == PROTRAITSisThereDict[ thisOrg ]:
#         concordance_score += 1 / NUM_ORGS

# # Are all the ones growing by KOMODO also growing by ProTraits? If yes, they may have a larger database and that may explain it.
# KOMODO_concordance_score = 0
# for thisOrg in prokNames:
#     if KOMODOisThereDict[ thisOrg ] or PROTRAITSisThereDict[ thisOrg ]:
#         KOMODO_concordance_score += 1

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
rxnDict = {}
for thisOrg in FULLPROTRAITSisThereDict.keys():
    rxnDict[ thisOrg ] = list( np.genfromtxt( 'organism_reactions/' + thisOrg + '.txt' ) )

# Comparing sizes.
sizes_ind = [ len( rxnDict[ thisOrgName ] ) 
              for thisOrgName in rxnDict.keys() 
              if FULLPROTRAITSisThereDict[ thisOrgName ] ] 
sizes_dep = [ len( rxnDict[ thisOrgName ] ) 
              for thisOrgName in rxnDict.keys() 
              if not FULLPROTRAITSisThereDict[ thisOrgName ] ] 
NUM_BINS = 15
ind_weights, dep_weights = [ 1 / len( sizes_ind ) ], [ 1 / len( sizes_dep ) ]
plt.hist( sizes_ind, bins=NUM_BINS, color='mediumseagreen', 
          histtype='stepfilled', normed=ind_weights )
plt.hist( sizes_dep, bins=NUM_BINS, color='dodgerblue', 
          histtype='stepfilled', normed=dep_weights )
plt.savefig( 'dists/sizes_ind_dep.svg' )
plt.show()

# Getting gene loss and gene gain lists.
ind_org_names = [ thisOrg for thisOrg in FULLPROTRAITSisThereDict.keys() 
                  if FULLPROTRAITSisThereDict[ thisOrg] ]
dep_org_names = [ thisOrg for thisOrg in FULLPROTRAITSisThereDict.keys() 
                  if not FULLPROTRAITSisThereDict[ thisOrg] ]
numGainsList, numLossesList = [], []
gainsList, lossesList = [], []

for thisIndOrg in tqdm( ind_org_names ):
    for thisDepOrg in tqdm( dep_org_names ):
        commonRxns = set( rxnDict[ thisIndOrg ] ) & set( rxnDict[ thisDepOrg ] )
        lostRxns = set( rxnDict[ thisIndOrg ] ).difference( commonRxns )
        gainRxns = set( rxnDict[ thisDepOrg ] ).difference( commonRxns )

        # Updating lists.
        gainsList.append( gainRxns )
        lossesList.append( lostRxns )
        numGainsList.append( len( gainRxns ) )
        numLossesList.append( len( lostRxns ) )

# Results: from ProTraits: mean number of gains = 304; losses = 591; mean gain/loss ratio = 0.5

# Plotting a gain/loss ratio distribution.
gainLossRatioList = [ numGainsList[ idx ] / numLossesList[ idx ] 
                      for idx in range( len( numGainsList ) ) 
                      if numLossesList[ idx ] and numGainsList[ idx ] ]f
ig, ax = plt.subplots(1)
bins= np.logspace( -4.5, 2, 20 )
plt.hist( gainLossRatioList, bins=bins, normed=True, histtype='stepfilled' )
ax.set_xscale('log')
plt.savefig( 'gain_to_loss_ratio_dist.svg' )
plt.show()

# Plotting a gain/loss ratio distribution.
plt.hist( numGainsList, bins=15, normed=True )
plt.hist( numLossesList, bins=15, normed=True )
plt.show()

# Getting gene loss and gene gain lists.
ind_org_names = [ thisOrg for thisOrg in FULLPROTRAITSisThereDict.keys() 
                  if FULLPROTRAITSisThereDict[ thisOrg ] ]
dep_org_names = [ thisOrg for thisOrg in FULLPROTRAITSisThereDict.keys() 
                  if not FULLPROTRAITSisThereDict[ thisOrg ] ]
numGainsList, numLossesList = [], []
gainsList, lossesList = [], []

for thisIndOrg in tqdm( FULLPROTRAITSisThereDict.keys() ):
    for thisDepOrg in tqdm( FULLPROTRAITSisThereDict.keys() ):
        if thisIndOrg != thisDepOrg:
            commonRxns = set( rxnDict[ thisIndOrg ] ) & set( rxnDict[ thisDepOrg ] )
            lostRxns = set( rxnDict[ thisIndOrg ] ).difference( commonRxns )
            gainRxns = set( rxnDict[ thisDepOrg ] ).difference( commonRxns )

            # Updating lists.
            gainsList.append( gainRxns )
            lossesList.append( lostRxns )
            numGainsList.append( len( gainRxns ) )
            numLossesList.append( len( lostRxns ) )

# Plotting a gain/loss ratio distribution.
            gainsList.append( gainRxns )
            lossesList.append( lostRxns )
            numGainsList.append( len( gainRxns ) )
gainLossRatioList = [ numGainsList[ idx ] / numLossesList[ idx ] 
                      for idx in range( len( numGainsList ) ) 
                      if numLossesList[ idx ] and numGainsList[ idx ] ]
ig, ax = plt.subplots(1)
bins= np.logspace( -4.5, 2, 20 )
plt.hist( gainLossRatioList, bins=bins, normed=True, histtype='stepfilled' )
ax.set_xscale('log')
plt.savefig( 'gain_to_loss_ratio_dist_any_to_any.svg' )
plt.show()

# Plotting a gain/loss ratio distribution.
plt.hist( numGainsList, bins=15, normed=True )
plt.hist( numLossesList, bins=15, normed=True )
plt.show()

# Checking for the tag 'free-living' in ProTraits and how it matches.
invOrgNameDict = { value: key for key, value in orgNameDict.items() }
isFreeLivingDict = {}
proTraits_freeliving_df = pd.read_csv( 'free_living_protraits.csv' )
for idx, thisRow in proTraits_freeliving_df.iterrows():
    try:
        thisOrgName = invOrgNameDict[ ' '.join( thisRow[0].split( )[ : 2 ] ) ]
    except:
        continue

    if bool( thisRow[-1] == 1 ):
        isFreeLivingDict[ thisOrgName ] = 1
    else:
        try:
            isFreeLivingDict[ thisOrgName ]
        except:
            isFreeLivingDict[ thisOrgName ] = 0

freeliving_concordance_score = 0
for thisOrg in FULLPROTRAITSisThereDict.keys():
    if isFreeLivingDict[ thisOrg ] and FULLPROTRAITSisThereDict[ thisOrg ]:
        print( orgNameDict[ thisOrg ] )
        freeliving_concordance_score += 1
