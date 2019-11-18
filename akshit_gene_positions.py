from ete3 import Tree
import numpy as np
import pandas as pd
from uniqify import uniqify
from unlistify import unlistify
from copy import copy
import pickle
import re
from load_kegg import *
from give_first_reactions import *

allByps = pd.read_csv('secreted_mets_segre.csv')
allNuts = pd.read_csv('seeds_from_vitkup.csv')
nuts = set(allNuts.iloc[:, 1].values)
byps = set(allByps.iloc[:, 1].values)
seedVec = np.zeros(len(rxnMat.T))
media = list(nuts | byps)
media = list(set(media) & set(cpds))
seedVec[[kegg_to_id[e] for e in Currency + media]] = 1

def firstReactions(genotype):
    kos = [pos_to_gene_map[i] for i in np.where((genotype > 0.5) * 1)[0]]
    refDF = pd.read_csv('KOREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])
    def give_number( string ):
        return int( string[-5:] )
    refDF['id'] = refDF['id'].apply( give_number )
    refDF['rid'] = refDF['rid'].apply( give_number )
    kos_to_rxns = refDF.loc[refDF['id'].isin(kos),:]
    can = refDF.loc[refDF['id'].isin(kos),:]['rid'].values
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in can]
    can = list((set(can)) & set(rxns))
    avrxns = np.array([rxn_kegg_to_id[e] for e in can])
    avScopeRxn = giveScope(rxnMat[avrxns], prodMat[avrxns], 
                           seedVec, sumRxnVec[avrxns])[1]
    scopeRxns = np.nonzero(avrxns * avScopeRxn)[0]
    starts = [kegg_to_id[c] for c in media]
    nuts_sum = np.sum(rxnMat[scopeRxns][:, starts], axis=1)
    firsts = [rxn_id_to_kegg[e] for e in avrxns[np.where(nuts_sum)[0]]]
    mfirsts = [int(l[1:]) for l in firsts]
    return uniqify(kos_to_rxns.loc[kos_to_rxns['rid'].isin(mfirsts)]['id'].tolist())

def lastReactions(genotype):
    kos = [pos_to_gene_map[i] for i in np.where((genotype > 0.5) * 1)[0]]
    refDF = pd.read_csv('KOREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])
    def give_number( string ):
        return int( string[-5:] )
    refDF['id'] = refDF['id'].apply( give_number )
    refDF['rid'] = refDF['rid'].apply( give_number )
    kos_to_rxns = refDF.loc[refDF['id'].isin(kos),:]
    can = refDF.loc[refDF['id'].isin(kos),:]['rid'].values
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in can]
    can = list((set(can)) & set(rxns))
    avrxns = np.array([rxn_kegg_to_id[e] for e in can])
    avrxns = np.array([rxn_kegg_to_id[e] for e in can])
    avScopeRxn = giveScope(rxnMat[avrxns], prodMat[avrxns], 
                           seedVec, sumRxnVec[avrxns])[1]
    scopeRxns = np.nonzero(avrxns * avScopeRxn)[0]
    cores = [kegg_to_id[c] for c in Core if c not in Currency]
    core_sum = np.sum(prodMat[scopeRxns][:, cores], axis=1)
    lasts = [rxn_id_to_kegg[e] for e in avrxns[np.where(core_sum)[0]]]
    mlasts = [int(l[1:]) for l in lasts]
    return uniqify(kos_to_rxns.loc[kos_to_rxns['rid'].isin(mlasts)]['id'].tolist())

def intReactions(genotype, firsts, lasts):
    kos = [pos_to_gene_map[i] for i in np.where((genotype > 0.5) * 1)[0]]
    refDF = pd.read_csv('KOREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])
    def give_number( string ):
        return int( string[-5:] )
    refDF['id'] = refDF['id'].apply( give_number )
    refDF['rid'] = refDF['rid'].apply( give_number )
    kos_to_rxns = refDF.loc[refDF['id'].isin(kos),:]
    can = refDF.loc[refDF['id'].isin(kos),:]['rid'].values
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in can]
    can = list((set(can)) & set(rxns))
    avrxns = np.array([rxn_kegg_to_id[e] for e in can])
    avScopeRxn = giveScope(rxnMat[avrxns], prodMat[avrxns], 
                           seedVec, sumRxnVec[avrxns])[1]
    scopeRxns = np.nonzero(avrxns * avScopeRxn)[0]
    ints = [rxn_id_to_kegg[e]  for e in avrxns[scopeRxns]]
    mints = [int(l[1:]) for l in ints]
    nkos = uniqify(kos_to_rxns.loc[kos_to_rxns['rid'].isin(mints)]['id'].tolist())
    kos_to_rxns = refDF.loc[refDF['id'].isin(nkos),:]
    return [e for e in kos_to_rxns['id'].tolist() 
            if e not in list(firsts) + list(lasts)]

gain_positions, loss_positions = [], []
for thisNode in tqdm(nodes):
    for thisChild in thisNode.children:
        gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
        loss_firsts = set(firstReactions(thisNode.genotype))
        loss_lasts = set(lastReactions(thisNode.genotype))
        loss_ints = set(intReactions(thisNode.genotype, loss_firsts, loss_lasts))

        gain_firsts = set(firstReactions(thisChild.genotype))
        gain_lasts = set(lastReactions(thisChild.genotype))
        gain_ints = set(intReactions(thisChild.genotype, gain_firsts, gain_lasts))

        filt_loss_firsts = loss_firsts - (loss_lasts | loss_ints)
        filt_loss_lasts = loss_lasts - (loss_firsts | loss_ints)
        filt_loss_ints = loss_ints - (loss_firsts | loss_lasts)

        filt_gain_firsts = gain_firsts - (gain_lasts | gain_ints)
        filt_gain_lasts = gain_lasts - (gain_firsts | gain_ints)
        filt_gain_ints = gain_ints - (gain_firsts | gain_lasts)

        gain_positions.append([])
        for this_gene in gainGenes:
            # if this_gene in ex_gains:
                if this_gene in filt_gain_firsts:
                    gain_positions[-1].append(0)
                elif this_gene in filt_gain_lasts:
                    gain_positions[-1].append(2)
                elif this_gene in filt_gain_ints:
                    gain_positions[-1].append(1)

        # if not gain_positions[-1]:
        #     gain_positions = gain_positions[:-1]
        # else:
        #     gain_positions[-1] = max(Counter(gain_positions[-1]), key=Counter(gain_positions[-1]).get)

        loss_positions.append([])
        for this_gene in lostGenes:
            # if this_gene in ex_losses:
                if this_gene in filt_loss_firsts:
                    loss_positions[-1].append(0)
                elif this_gene in filt_loss_lasts:
                    loss_positions[-1].append(2)
                elif this_gene in filt_loss_ints:
                    loss_positions[-1].append(1)

        # if not loss_positions[-1]:
        #     loss_positions = loss_positions[:-1]
        # else:
        #     loss_positions[-1] = max(Counter(loss_positions[-1]), key=Counter(loss_positions[-1]).get)


from collections import Counter
gain_tallies = [Counter(gain_positions[i]) for i in range(len(gain_positions)) 
                 if gain_positions[i]]
loss_tallies = [Counter(loss_positions[i]) for i in range(len(loss_positions))
                if loss_positions[i]]


frac_gain_tallies = [{g : (event[g] / sum(event.values())) for g in event} 
                     for event in gain_tallies]
frac_loss_tallies = [{g : (event[g] / sum(event.values())) for g in event} 
                     for event in loss_tallies]

def meanAtPos(frac_tallies, pos):
    all_fracs = []
    for i in range(len(frac_tallies)):
        try:
            all_fracs.append(frac_tallies[i][pos])
        except:
            pass
            all_fracs.append(0)
    return np.mean(all_fracs)

def medianAtPos(frac_tallies, pos):
    all_fracs = []
    for i in range(len(frac_tallies)):
        try:
            all_fracs.append(frac_tallies[i][pos])
        except:
            pass
            all_fracs.append(0)
    return np.median(all_fracs)

mean_gain_tallies = {pos: meanAtPos(frac_gain_tallies, pos)
                    for pos in [0, 1, 2]}
mean_loss_tallies = {pos: meanAtPos(frac_loss_tallies, pos)
                    for pos in [0, 1, 2]}

median_gain_tallies = {pos: medianAtPos(frac_gain_tallies, pos)
                    for pos in [0, 1, 2]}
