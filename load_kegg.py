import numpy as np
import pickle

#-------------------------------------------------------------------------
# Loading all relevant data files.
#-------------------------------------------------------------------------
stoich_matrix = np.load('../pangenome_cross_feeding/akshit_stoich_matrix.npy')
names = [e[:-1] for e in open('../pangenome_cross_feeding/akshit_names.txt', 'r').readlines()]
rxns = ''.join(open('../pangenome_cross_feeding/akshit_rxns.txt', 'r').readlines()).split()
cpds = ''.join(open('../pangenome_cross_feeding/akshit_mets.txt', 'r').readlines()).split()

met_map = kegg_to_id = {cpds[e] : e for e in range(len(cpds))}
inv_met_map = id_to_kegg = {value: key for key, value in met_map.items()}
rxn_map = rxn_kegg_to_id = {rxns[e] : e for e in range(len(rxns))}
inv_rxn_map = rxn_id_to_kegg = {value: key for key, value in rxn_map.items()}
cpd_string_dict = {cpds[e] : names[e] for e in range(len(names))}

#-------------------------------------------------------------------------
# Initializing and setting things up.
#-------------------------------------------------------------------------
Currency = ''.join(open('../pangenome_cross_feeding/bigg_currency.txt', 'r').readlines()).split()
Core = ''.join(open('../pangenome_cross_feeding/bigg_core.txt', 'r').readlines()).split()

#-------------------------------------------------------------------------
# Defining reactant, product and reactant number vectors for the scope 
# expansion algorithm.
#-------------------------------------------------------------------------
rxnMat = ((stoich_matrix.clip(max = 0.0)) != 0) * 1
prodMat = (stoich_matrix.clip(min = 0.0) != 0) * 1
sumRxnVec = np.sum(rxnMat, axis = 1)
sumProdVec = np.sum(prodMat, axis = 1)
