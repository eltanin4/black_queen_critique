import numpy as np
import pickle

NUM_ORGS = 835
# Getting the list of important 'terminal' metabolites.
mod_to_comp_dict = pickle.load(open('module_containers/mod_to_comp_dict.dat', 'rb'))
terMets = list(set([ mod_to_comp_dict[thisMod][-1] for thisMod in mod_to_comp_dict 
                     if mod_to_comp_dict[thisMod] ]))
terMets.remove('CHPF2;')  # Spurious metabolite being removed. What the fuck is this even?

#-------------------------------------------------------------------------
# Loading all relevant data files.
#-------------------------------------------------------------------------
stoich_matrix = np.genfromtxt('segre_stoich_matrix.csv', delimiter=',').T
names = [e[:-1] for e in open('segre_names.txt', 'r').readlines()]
cpds = [e[:-1] for e in open('segre_mets.txt', 'r').readlines()]
rxns = [e[:-1] for e in open('segre_rxns.txt', 'r').readlines()]

met_map = kegg_to_id = {cpds[e] : e for e in range(len(cpds))}
inv_met_map = id_to_kegg = {value: key for key, value in met_map.items()}
rxn_map = rxn_kegg_to_id = {rxns[e] : e for e in range(len(rxns))}
inv_rxn_map = rxn_id_to_kegg = {value: key for key, value in rxn_map.items()}
indices_of_currency_mets = np.genfromtxt(
    'kegg_curr.csv', delimiter=',')
cpd_string_dict = {cpds[e] : names[e] for e in range(len(names))}

#-------------------------------------------------------------------------
# Initializing and setting things up.
#-------------------------------------------------------------------------
Core = list(set(terMets[:]) & set(cpds))
Currency = ['C' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in indices_of_currency_mets]

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
