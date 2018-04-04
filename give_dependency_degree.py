import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify


NUM_ORGS = 835


# Getting the list of important 'terminal' metabolites.
mod_to_comp_dict = pickle.load(open('module_containers/mod_to_comp_dict.dat', 'rb'))
terMets = list(set([ mod_to_comp_dict[thisMod][-1] for thisMod in mod_to_comp_dict 
                     if mod_to_comp_dict[thisMod] ]))
terMets.remove('CHPF2;')  # Spurious metabolite being removed. What the fuck is this even?


# Note that glycans have been ignored in this analysis. 
# We get a set of 228 terminal metabolites thereby.


# Getting the set of bacterial abbreviations for each organism in KEGG.
orgNames = []
with open('endo_removed_prok_abbr_kegg.txt', 'r') as f:
    for thisLine in f.readlines():
        orgNames.append( thisLine.strip() )


# Getting the compound to reaction mapping.
cpd_to_rxn_df = pd.read_csv('compound_to_reaction.txt', 
                            delimiter='\t', sep='delimiter', 
                            header=None, names=['cpd', 'rxn'])


# Now calculating a 'producibility' score for each organism.
dep_degree = {}
for thisOrg in tqdm(orgNames):

  # Getting all the reactions performable by this organism.
  dep_degree[thisOrg] = {}
  can = set(np.genfromtxt('organism_reactions/' + thisOrg + '.txt').astype(int))

  for thisMet in tqdm(terMets):
    # Need all reactions that this metabolite is part of.
    relRxns = list(cpd_to_rxn_df.loc[
                   cpd_to_rxn_df['cpd'] == 'cpd:' + thisMet]['rxn'].values)

    # Converting this to a set of ints.
    allowed = set([int(e) for e in re.findall(r'\d+', ''.join(relRxns))])

    # Saving whether this organism can make this metabolite.
    dep_degree[thisOrg][thisMet] = bool(allowed & can)

# Saved dep_degree and dep_degree_dict.dat



# Now plotting a distribution of producibility scores across organisms.
indep_scores = {thisOrg : sum(dep_degree[thisOrg].values()) for thisOrg in orgNames}
num_genes = {thisOrg : len(set(np.genfromtxt('organism_kogenes/' + thisOrg + '.txt').astype(int)))
             for thisOrg in orgNames}

genes, indscore = [], []
for thisOrg in orgNames:
  genes.append(num_genes[thisOrg]) 
  indscore.append(indep_scores[thisOrg])
fig, ax = plt.subplots(1)
ax.scatter(genes, indscore)
plt.show()
