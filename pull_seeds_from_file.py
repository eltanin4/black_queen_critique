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
