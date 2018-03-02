import numpy as np
import pickle

NUM_ORGS = 435

 # Getting all organism seed strings
orgSeedList = []
 with open('seeds.txt', 'r') as f:
     for thisOrg in range(NUM_ORGS):              
         thisLine = f.readline()    
         orgSeedList.append(thisLine)

def find_between_r(s, first, last):
    try:
        start = s.rindex(first) + len(first)
        end = s.rindex(last, start)
        return s[start:end]
    except ValueError:
        return ""

# Converting all strings to list of list of numbers, each ID one of the KEGG seed IDs.
orgSeedList = [ [ int(e) for e in find_between_r( orgSeed, ':', '\n').split(' ')[1:-1] ] for orgSeed in orgSeedList ]

pickle.dump( orgSeedList, open( 'org_seed_list.dat', 'wb' ) )

# Getting sizes for each organism.
rxnDict = {}
for thisOrg in prokNames:
    rxnDict[ thisOrg ] = list( np.genfromtxt( 'organism_reactions/' + thisOrg + '.txt' ) )

# Get the dependency score distribution and plotting it with some standard organism values.
depScores = list( pickle.load( open( 'dict_org_dep_scores.dat', 'rb' ) ).values() )
normedDepWeights = np.ones_like( depScores ) / len( depScores )
fig, ax = plt.subplots(1)
plt.hist( depScores, bins=40, histtype='stepfilled', weights=normedDepWeights )
ax.axvline( depScoreDict[ 'eco' ], color='black', linewidth=2.5)
plt.savefig( 'dists/init_dep_score_dict_with_eco.svg' )
plt.show()

# Getting the names of the abbrevations from KEGG for KOMODO comparison.
orgNameDict = {}
with open('prok_kegg_abbr_with_names', 'r') as f:
    for thisOrg in range(NUM_ORGS):              
        thisLine = f.readline()    
        orgNameDict[ thisLine[:3] ] = thisLine[4:-1]

prokNames = list( orgNameDict.keys() ) 

# March 1, 2018
MEM_CUTOFF = 0.95
seedMembershipDict = {}
all_seeds = uniqify( unlistify( orgSeedList ) )
for thisSeed in all_seeds:
    seedMembershipDict[ thisSeed ] = sum( [ bool( thisSeed in orgSeedDict[ thisOrg ] ) 
                                          for thisOrg in prokNames ] )

# Seed membership histogram is bimodal (U-shaped), so extracting top 10% of seed compounds.
PERCENT_CHOSEN = 10
listSeedMembership = list( seedMembershipDict.values() )
coreSeeds = sorted( range( len( listSeedMembership ) ), 
                    key=lambda i: listSeedMembership[ i ], reverse=True )[ : 
                    int( PERCENT_CHOSEN / 100 * len( prokNames ) ) ]

# Now removing all core seeds to make a non-core seed dict for all organisms.
orgSeedDict_CORE = {}
for thisOrg in prokNames:
    orgSeedDict_CORE[ thisOrg ] = [ thisSeed for thisSeed in orgSeedDict[ thisOrg ] 
                                    if thisSeed not in coreSeeds]

# Now getting the "core seed normalized" dependency scores.
depScoreDict_CORE = {}
for thisOrg in prokNames:
    depScoreDict_CORE[ thisOrg ] = ( len( orgSeedDict_CORE[ thisOrg ] ) / 
                                     len( orgRxnDict[ thisOrg ] ) )

# Now plotting all core normalized dependency scores.
depScores_CORE = list( depScoreDict_CORE.values() )
depScores_COREWeights = np.ones_like( np.array( depScores_CORE ) ) / len( depScores_CORE )
plt.hist( depScores_CORE, histtype='stepfilled', bins=40, weights=depScores_COREWeights )
plt.show()

# Getting ProTraits grower list.
proGrowersList = []
with open('proTraits_growers.txt', 'r') as f:
    for thisLine in f:
        proGrowersList.append( thisLine[:-1] )

# Checking against pro-traits growers on many many media.
isThereDict = {}
for thisOrg in prokNames:
    isThereDict[ thisOrg ] = bool( orgNameDict[ thisOrg ] in proGrowersList )

# Checking dependency scores for all ProTraits predicted organisms.
indDepScoreDict = { thisOrg: depScoreDict[ thisOrg ] 
                    for thisOrg in prokNames 
                    if  isThereDict[ thisOrg ] }

# March 2, 2018
# Getting KOMODO growers again.
import pandas as pd
org_to_media_df = 

