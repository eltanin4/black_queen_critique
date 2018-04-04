import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify
from load_kegg import *
from give_scope import giveScope

# Getting the set of bacterial abbreviations for each organism in KEGG.
orgNames = []
with open('endo_removed_prok_abbr_kegg.txt', 'r') as f:
    for thisLine in f.readlines():
        orgNames.append( thisLine.strip() )
        

# Generating a vector of all metabolites initially provided, i.e. seeds.
seedVec = np.zeros(len(rxnMat.T))
seeds_df = pd.read_csv('seeds_from_vitkup.csv')
media = list(seeds_df['kegg_id'].values)
media = list(set(media) & set(cpds))
seedVec[[kegg_to_id[e] for e in Currency + media]] = 1

# Now calculating a 'producibility' score for each organism.
ind_degree = {}
for thisOrg in tqdm(orgNames):
    # Getting all the reactions performable by this organism.
    can = set(np.genfromtxt('organism_reactions/' + thisOrg + '.txt').astype(int))
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in can]
    can = list((set(can)) & set(rxns))
    avrxns = [rxn_kegg_to_id[e] for e in can]

    # Calculating the metabolites within the scope of 
    # this organism's reaction network. 
    scopeMets = giveScope(rxnMat[avrxns], prodMat[avrxns], seedVec, sumRxnVec[avrxns])[0]
    
    # Finding how much of the core is within the network's scope.
    ind_degree[thisOrg] = list(set([id_to_kegg[e] for e in set(np.where(scopeMets)[0])]) & set(Core))

# Saved ind_degree and ind_degree_dict.dat


# Now plotting a distribution of producibility scores across organisms.
indep_scores = {thisOrg : len(ind_degree[thisOrg]) for thisOrg in orgNames}
num_genes = {thisOrg : len(set(np.genfromtxt('organism_kogenes/' + thisOrg + '.txt').astype(int)))
             for thisOrg in orgNames}
num_rxns = {thisOrg : len(set(np.genfromtxt('organism_reactions/' + thisOrg + '.txt').astype(int)))
            for thisOrg in orgNames}

# Getting just a list of corresponding gene numbers, reactions numbers and ind_scores
genes, indscore, reactions = [], [], []
for thisOrg in orgNames:
  genes.append(num_genes[thisOrg])
  reactions.append(num_rxns[thisOrg])
  indscore.append(indep_scores[thisOrg])

# Normalizing the independence scores between 0 and 1.
indscore = np.asarray(indscore)
norm_indscore = (indscore - min(indscore)) / (max(indscore) - min(indscore))
