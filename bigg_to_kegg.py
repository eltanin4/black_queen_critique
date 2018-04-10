import pandas as pd
import numpy as np
from uniqify import uniqify
import urllib3
import urllib.request
import shutil
from tqdm import tqdm

handler = urllib.request.ProxyHandler({'http': 'proxy.ncbs.res.in:3128'})
opener = urllib.request.build_opener(handler)
urllib.request.install_opener(opener)
urllib.request.ProxyHandler({'ftp': 'proxy.ncbs.res.in:3128'})

#-------------------------------------------------------------------------
# Getting a filtered list of BiGG database metabolites. 
#-------------------------------------------------------------------------
bigg_df = pd.read_table('bigg_models_metabolites.csv')
hasDatabase = np.logical_not(pd.isnull(bigg_df['database_links']).values)
bigg_df = bigg_df.loc[np.where(hasDatabase)[0], :]


#-------------------------------------------------------------------------
# Getting all organism models in BiGG and identifying their biomass 
# reaction strings.
#-------------------------------------------------------------------------
# all_bigg_models = open('all_bigg_models.txt', 'r').readline()
# org_biggids = re.findall(r'\"bigg_id\": \"(.*?)\"', all_bigg_models)

# # Getting all organism reactions.
# saveDir = 'bigg_orgs/'
# reaction_url_holder = 'http://bigg.ucsd.edu/api/v2/models/'
# for thisOrg in tqdm(org_biggids):
#     thisURL = reaction_url_holder + thisOrg + '/reactions'
#     urllib.request.urlretrieve( thisURL, saveDir + thisOrg + '.txt' )

# Tracking the biomass reaction name in each organism's reaction set.
saveDir = 'bigg_orgs/'
biomass_url_holder = 'http://bigg.ucsd.edu/api/v2/models/'
org_bm_list = []
for thisOrg in tqdm(org_biggids):
    with open(saveDir + thisOrg + '.txt', 'r') as f:
        s = f.readline()
        thisBMname = re.findall(r'\"bigg_id\": \"(BIOMASS.*?)\"', s)[0]
        org_bm_list.append(thisBMname)
        # thisURL = biomass_url_holder + thisOrg + '/reactions/' + thisBMname
        # urllib.request.urlretrieve(thisURL, saveDir + thisBMname + '.txt')

# Getting the unique list.
org_bm_list = uniqify(org_bm_list)


#-------------------------------------------------------------------------
# Now mapping the metabolite names to possible KEGG IDs.
#-------------------------------------------------------------------------
bmKEGGS = []
for thisBMname in tqdm(org_bm_list):
    thisbm = open(saveDir + thisBMname + '.txt', 'r').readline()
    biggids = re.findall(r'\"bigg_id\": \"(.*?)\"', thisbm)[:-2]
    bmMetsDF = pd.DataFrame({'universal_bigg_id':biggids})
    filt_bigg_df = bigg_df.merge(bmMetsDF, on='universal_bigg_id')
    relLinks = list(filt_bigg_df['database_links'].values)

    # Now identifying all the KEGG IDs.
    for tl in relLinks:
        try:
            # Finding the relevant KEGG compound strings.
            bmKEGGS.append(re.findall(r'KEGG Compound: http://identifiers.org/kegg.compound/(.*?);', tl)[0])

            # Explicitly removing glycans.
            if bmKEGGS[-1][0] == 'G':
                bmKEGGS = bmKEGGS[:-1]
        except:
            pass

# Removing duplicates.
bmKEGGS = uniqify(bmKEGGS)
