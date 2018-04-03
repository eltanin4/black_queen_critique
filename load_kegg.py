import numpy as np
import pickle

#-------------------------------------------------------------------------
# Loading all relevant data files.
#-------------------------------------------------------------------------
stoich_matrix = np.genfromtxt('segre_stoich_matrix.csv', delimiter=',')
# stoich_matrix = np.genfromtxt(
#     '../new_kegg/data_files/random_stoich_matrix.txt')
met_map = kegg_to_id = pickle.load(
    (open('../new_kegg/data_files/pure_kegg_to_self_met_map.dat', 'rb')))
inv_met_map = id_to_kegg = {value: key for key, value in met_map.items()}
rxn_map = rxn_kegg_to_id = pickle.load(
    (open('../new_kegg/data_files/pure_kegg_to_self_rxn_map.dat', 'rb')))
inv_rxn_map = rxn_id_to_kegg = {value: key for key, value in rxn_map.items()}
indices_of_core_mets = np.genfromtxt(
    '../new_kegg/data_files/kegg_core.csv', delimiter=',')
# indices_of_core_mets = np.genfromtxt(
#       '../new_kegg/data_files/rand_core.csv', delimiter=',')
indices_of_currency_mets = np.genfromtxt(
    '../new_kegg/data_files/kegg_curr.csv', delimiter=',')
indices_of_energy_mets = np.genfromtxt(
    '../new_kegg/data_files/kegg_atp.csv', delimiter=',')
mets = list(met_map.values())
display_lookup = pickle.load(
    open('../new_kegg/data_files/rxn_string_dict.dat', 'rb'))
cpd_string_dict = pickle.load(
    open('../new_kegg/data_files/cpd_string_dict.dat', 'rb'))

#-------------------------------------------------------------------------
# Initializing and setting things up.
#-------------------------------------------------------------------------
Core = [kegg_to_id[i] for i in indices_of_core_mets]
Currency = [kegg_to_id[i] for i in indices_of_currency_mets]
Energy = [kegg_to_id[i] for i in indices_of_energy_mets]

#-------------------------------------------------------------------------
# Defining reactant, product and reactant number vectors for the scope 
# expansion algorithm.
#-------------------------------------------------------------------------
rxnMat = ((stoich_matrix.clip(max = 0.0)) != 0) * 1
prodMat = (stoich_matrix.clip(min = 0.0) != 0) * 1
sumRxnVec = np.sum(rxnMat, axis = 1)
sumProdVec = np.sum(prodMat, axis = 1)

#-------------------------------------------------------------------------
# Defining the clipped reactant-only and product-only matrices.
#-------------------------------------------------------------------------
rho = stoich_matrix.clip(max = 0.0)
pi = stoich_matrix.clip(min = 0.0)

# #-------------------------------------------------------------------------
# # Generating the anaerobic stoichiometric matrix.
# #-------------------------------------------------------------------------
# oxInvolvingRxns = np.where(stoich_matrix[:, kegg_to_id[7]] != 0.0)
# stoich_matrix[oxInvolvingRxns] = np.zeros(np.shape(stoich_matrix[oxInvolvingRxns])[1])

# #-------------------------------------------------------------------------
# # Defining reactant, product and reactant number vectors for the scope 
# # expansion algorithm.
# #-------------------------------------------------------------------------
# rxnMat = ((stoich_matrix.clip(max = 0.0)) != 0) * 1
# prodMat = (stoich_matrix.clip(min = 0.0) != 0) * 1
# sumRxnVec = np.sum(rxnMat, axis = 1)
# sumProdVec = np.sum(prodMat, axis = 1)

# Currency = [met for met in Currency if met != kegg_to_id[7]]
