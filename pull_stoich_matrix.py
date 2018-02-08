import os
import numpy as np
from bs4 import BeautifulSoup
import urllib3
import urllib.request
import shutil
import pickle
import re
import collections
from print_progress_bar import *

rxnids = np.genfromtxt('all_reaction_ids_kegg.txt', delimiter='\n')
stoich_matrix = np.zeros((len(rxnids) * 2, 99999))

saveDir = 'reaction_files/'
if not os.path.exists( saveDir ):
    os.mkdir( saveDir )

#-------------------------------------------------------------------------
# This code extracts all KEGG webpages in the pathway map 01120.
#-------------------------------------------------------------------------
# c = urllib3.PoolManager()
rxn_url_holder = 'http://www.genome.jp/dbget-bin/www_bget?rn:R'

handler = urllib.request.ProxyHandler({'http': 'proxy.ncbs.res.in:3128'})
opener = urllib.request.build_opener(handler)
urllib.request.install_opener(opener)
urllib.request.ProxyHandler({'ftp': 'proxy.ncbs.res.in:3128'})

for i in range(len(rxnids)):
    print_progress_bar(i, len(rxnids), 'Pulling all reaction files')
    r_id = str(int(rxnids[i]))
    current_url = rxn_url_holder + '0' * (5 - len(r_id)) + r_id
    path = 'reaction_files/' + '0' * (5 - len(r_id)) + r_id + '.html'
    if not os.path.isfile(path):
        urllib.request.urlretrieve(current_url, path)


def find_between_r(s, first, last):
    try:
        start = s.rindex(first) + len(first)
        end = s.rindex(last, start)
        return s[start:end]
    except ValueError:
        return ""

exceptions = [2362, 2887, 5196]
rxn_string_dict = {}
#-------------------------------------------------------------------------
# This converts the extracted html files into reactions and then into a bloated
# stoichiometric matrix.
#-------------------------------------------------------------------------
irreversibles = []
for i in range(len(rxnids)):
    print_progress_bar(i, len(rxnids), 'Pulling reactions from files')
    r_id = str(int(rxnids[i]))
    current_url = '0' * (5 - len(r_id)) + r_id
    path = 'reaction_files/' + '0' * (5 - len(r_id)) + r_id + '.html'

    soup = BeautifulSoup(open(path, 'rb'))
    letters = soup.find_all('td', class_='td21')
    letters_2 = soup.find_all('td', class_='td20')
    for el in letters + letters_2:
        if 'cpd' in str(el):
            eq = str(el)
            break

    a = find_between_r(str(soup), 'Definition',
                       '<br/>')[:1500].split('br/>')[0]
    rxn_string = find_between_r(a, 'n\">', '<')
    rxn_string_dict[rxnids[i]] = rxn_string

    # Checking if the reaction is reversible or irreversible.
    is_reversible = True
    try:
        reactant_side, product_side = re.split(r'\s*&lt;=&gt;\s*', eq)
    except:
        is_reversible = False
        reactant_side, product_side = re.split(r'\s*=&gt;\s*', eq)

    if not int(r_id) in exceptions:
        reactant_stoichs = re.findall(r'n\">(.*?) <a href', reactant_side)
        reactant_stoichs += re.findall(r'a> +(.*?) <a href', reactant_side)
        # reactant_stoichs = [el for el in reactant_stoichs if 'C' not in el]
        reactants = re.findall(r'a href=\"(.*?)\"', reactant_side)
        reactant_ids = [int(r[-5:]) for r in reactants]

        reactant_dict, flag = {}, 0
        if reactant_stoichs == []:
            reactant_dict[reactant_ids[0]] = 1

        elif '+' in reactant_stoichs[0]:
            reactant_dict[reactant_ids[0]] = 1
        else:
            reactant_dict[reactant_ids[0]] = int(reactant_stoichs[0])
            flag = 1
        for x, elem in enumerate(reactant_stoichs[1:]):
            if not flag:
                for q in elem.split(' '):
                    try:
                        num = int(q)
                    except:
                        num = -1
                if num != -1:
                    reactant_dict[reactant_ids[x + 1]] = num
                elif num == -1:
                    reactant_dict[reactant_ids[x + 1]] = 1
            elif flag:
                for q in elem.split(' '):
                    try:
                        num = int(q)
                    except:
                        num = -1
                if num != -1:
                    reactant_dict[reactant_ids[x + 1]] = num
                elif num == -1:
                    reactant_dict[reactant_ids[x + 1]] = 1
        if len(list(reactant_dict.keys())) != len(reactants):
            reactant_dict[reactant_ids[-1]] = 1
    else:
        reactant_dict = {some_id: 1 for some_id in reactant_ids}

    if not int(r_id) in exceptions:
        product_stoichs = re.findall(r'(.*?) <a href', product_side)
        product_stoichs += re.findall(r'a> +(.*?) <a href', product_side)
        product_stoichs = [el for el in product_stoichs if 'C' not in el]
        products = re.findall(r'a href=\"(.*?)\"', product_side)
        product_ids = [int(r[-5:]) for r in products]

        product_dict, flag = {}, 0
        if len(product_stoichs) != len(product_ids):
            product_stoichs = [''] + product_stoichs
        if '+' in product_stoichs[0]:
            product_dict[product_ids[0]] = 1
        elif product_stoichs[0] == '':
            product_dict[product_ids[0]] = 1
            # flag = 1
        else:
            product_dict[product_ids[0]] = int(product_stoichs[0])
            flag = 1
        for x, elem in enumerate(product_stoichs[1:]):
            if not flag:
                for q in elem.split(' '):
                    try:
                        num = int(q)
                    except:
                        num = -1
                if num != -1:
                    product_dict[product_ids[x + 1]] = num
                elif num == -1:
                    product_dict[product_ids[x + 1]] = 1
            elif flag:
                for q in elem.split(' '):
                    try:
                        num = int(q)
                    except:
                        num = -1
                if num != -1:
                    product_dict[product_ids[x + 1]] = num
                elif num == -1:
                    product_dict[product_ids[x + 1]] = 1
        if len(list(product_dict.keys())) != len(products):
            product_dict[product_ids[-1]] = 1
    else:
        product_dict = {some_id: 1 for some_id in product_ids}

    for this_id in reactant_ids:
        stoich_matrix[i][this_id] = -1 * reactant_dict[this_id]
    for this_id in product_ids:
        stoich_matrix[i][this_id] = +1 * product_dict[this_id]
    if is_reversible:
        for this_id in reactant_ids:
            stoich_matrix[len(rxnids) + i][this_id] = + \
                1 * reactant_dict[this_id]
        for this_id in product_ids:
            stoich_matrix[len(rxnids) + i][this_id] = - \
                1 * product_dict[this_id]
    else:
        irreversibles.append(i)

print('\n')

# Now I have a bloated stoichiometric matrix that needs to shedd off a lot
# of stuff.
transposed_matrix = stoich_matrix.transpose()
empty_mets = []
for met_id, met in enumerate(transposed_matrix):
    print_progress_bar(met_id, len(transposed_matrix),
                       'Figuring out unused metabolites')
    if met.any():
        pass
    else:
        empty_mets.append(met_id)

print('\n')

# Now for mapping KEGG IDs to my personal IDs so I can use this stuff.
used_met_ids = [m for m in list(range(99999)) if m not in empty_mets]
kegg_to_self_met_map = {}
new_transpose = []
for nt_counts, met_id in enumerate(used_met_ids):
    print_progress_bar(met_id, len(used_met_ids),
                       'Removing unused metabolites')
    kegg_to_self_met_map[met_id] = nt_counts
    new_transpose.append(transposed_matrix[met_id])
stoich_matrix = np.asarray(new_transpose).transpose()

kegg_to_self_rxn_map = {}
reduced_stoich_matrix = []
for r_id, rxn in enumerate(stoich_matrix):
    if rxn.any():
        reduced_stoich_matrix.append(rxn)
        if r_id < len(rxnids):
            kegg_to_self_rxn_map[int(rxnids[r_id])] = len(
                reduced_stoich_matrix) - 1
stoich_matrix = np.asarray(reduced_stoich_matrix)
print('\n')

# types = []
# for r_ind, reaction in enumerate(stoich_matrix):
#     print_progress_bar(r_ind, len(stoich_matrix),
#                        'Pilfering reactions in types')
#     counts = dict(collections.Counter(reaction))
#     try:
#         rtype = (num_reactants, num_products) = (counts[-1], counts[+1])
#         types.append(rtype)
#     except:
#         pass

# print('\n')

# rcounts = dict(collections.Counter(types))
# filtered_types = {i: rcounts[i] for i in rcounts.keys() if rcounts[i] <= 1000}
# allowed_types = list(filtered_types.keys())

# #-------------------------------------------------------------------------
# # Generating the reaction distribution.
# #-------------------------------------------------------------------------
# reduced_stoich_matrix = []
# for r_ind, reaction in enumerate(stoich_matrix):
#     print_progress_bar(r_ind, len(stoich_matrix),
#                        'Generating reaction distribution')
#     counts = dict(collections.Counter(reaction))
#     try:
#         rtype = (num_reactants, num_products) = (counts[-1], counts[+1])
#         if rtype in allowed_types:
#             reduced_stoich_matrix.append(reaction)
#     except:
#         pass
# stoich_matrix = np.asarray(reduced_stoich_matrix)

# reaction_density = np.zeros((max(np.asarray(allowed_types)[:, 0]) + 1,
#                              max(np.asarray(allowed_types)[:, 1]) + 1))
# for i, j in filtered_types:
#     reaction_density[i, j] = filtered_types[(i, j)]
# reaction_density /= np.sum(reaction_density)

# plot_heatmap(reaction_density)