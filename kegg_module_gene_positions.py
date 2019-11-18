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
import pandas as pd

def find_between_r(s, first, last):
    try:
        start = s.rindex(first) + len(first)
        end = s.rindex(last, start)
        return s[start:end]
    except ValueError:
        return ""

#-------------------------------------------------------------------------
# First getting a hierarchical list of all KEGG module categories 
# and names.
#-------------------------------------------------------------------------
category_dict, module_list = { }, [ ]
with open( 'kegg_pathway_module_htext.txt', 'r' ) as modFile:
    # Going over the htext line-by-line.
    for thisLine in modFile.readlines():
        # Checking if new category reached.
        if thisLine[0] == 'C':
            catName = find_between_r( thisLine, '  ', '\n' )
            category_dict[ catName ] = []

        # Now filling all modules in this category.
        if thisLine[0] == 'D':
            category_dict[ catName ].append( find_between_r( thisLine, '   ', '\n' ) )
            module_list.append( find_between_r( thisLine, '   ', '\n' ) )

#-------------------------------------------------------------------------
# This code extracts all KEGG webpages for KEGG module data.
#-------------------------------------------------------------------------
module_url_holder = 'http://www.genome.jp/kegg-bin/show_module?'

handler = urllib.request.ProxyHandler({'http': 'proxy.ncbs.res.in:3128'})
opener = urllib.request.build_opener(handler)
urllib.request.install_opener(opener)
urllib.request.ProxyHandler({'ftp': 'proxy.ncbs.res.in:3128'})

# Getting the overall save directory.
saveDir = 'module_exprs/'
if not os.path.exists( saveDir ):
    os.mkdir( saveDir )

# Getting the module names.
module_list = pickle.load( open( 'module_list.dat', 'rb' ) )

# Pulling all module pages in KEGG MODULE.
for modIndex, modName in enumerate( module_list ):
    print_progress_bar( modIndex, len( module_list ), 'Getting all module pages')
    current_url = module_url_holder + str( modName[ : 6 ] )
    path = saveDir + str( modName[ : 6 ] ) + '.html'
    if not os.path.isfile(path):
        urllib.request.urlretrieve( current_url, path )

# Getting the logical expression string corresponding to each module.
expr_dict = { }
for modIndex, modName in enumerate( module_list ):
    print_progress_bar( modIndex, len( module_list ), 'Getting all logical expression strings')
    # Get the list of strings, each row corresponding to an EC number.
    try:
        with open( saveDir + str( modName[ : 6 ] ) + '.html', 'r' ) as modFile:
            modString = BeautifulSoup(modFile, "lxml").text
    except:
        continue

    # Extracting the right expression string.
    extracted_str = find_between_r( modString, 'Definition', 'Type\n' )
    preproc_str = re.sub( '\s+', ' ', extracted_str ).strip()
    final_str = preproc_str.replace( '--', '' ).strip()
    expr_dict[ modName[ :6 ] ] = final_str

# The expr_dict is now saved in the module_exprs directory.

#-------------------------------------------------------------------------
# This gets all KO and compound IDs associated with the KEGG module.
#-------------------------------------------------------------------------
module_url_holder = 'http://www.genome.jp/kegg-bin/module_ko_list?map='

# Getting the overall save directory.
saveDir = 'module_containers/'
if not os.path.exists( saveDir ):
    os.mkdir( saveDir )

# Getting the module names.
module_list = pickle.load( open( 'module_list.dat', 'rb' ) )

# Pulling all module pages in KEGG MODULE.
for modIndex, modName in enumerate( module_list ):
    print_progress_bar( modIndex, len( module_list ), 'Getting all module pages')
    current_url = module_url_holder + str( modName[ : 6 ] ) + '&org=ko'
    path = saveDir + str( modName[ : 6 ] ) + '.html'
    if not os.path.isfile(path):
        urllib.request.urlretrieve( current_url, path )

# Getting the KO genes and compounds which are in the module. 
mod_to_gene_dict, mod_to_comp_dict = { }, { }
for modIndex, modName in enumerate( module_list ):
    print_progress_bar( modIndex, len( module_list ), 'Getting all KO genes and compounds')
    # Get the list of strings, each row corresponding to an EC number.
    try:
        with open( saveDir + str( modName[ : 6 ] ) + '.html', 'r' ) as modFile:
            modString = BeautifulSoup(modFile, "lxml").text
    except:
        continue

    # Extracting the right expression string.
    extracted_str = modString[:]
    preproc_str = re.sub( '\s+', ' ', extracted_str ).strip()
    geneNames = [ t.replace( '-' , '0' ) for t in preproc_str.split() if t.startswith('K') and len(t) == 6 ]
    compNames = [ t.replace( '-' , '0' ) for t in preproc_str.split() if t.startswith('C') and len(t) == 6 ]
    mod_to_gene_dict[ modName[ : 6 ] ] = geneNames
    mod_to_comp_dict[ modName[ : 6 ] ] = compNames

# Saving mod_to_gene_dict
# Note: the above mod_to_gene map is flawed. It's more reliable to extract it from downstairs.

# Clearing the expr_dict
import re
expr_dict = pickle.load( open( 'module_exprs/expr_dict.dat', 'rb' ) )
old_expr_dict = deepcopy( expr_dict )
for thisMod in expr_dict:
    expr_dict[ thisMod ] = re.sub( '(K\d+)(K\d+)', '\g<2>', expr_dict[ thisMod ] )
    expr_dict[ thisMod ] = expr_dict[ thisMod ].replace( ' K', 'K' ).strip().replace( '  ', ' ' )
    expr_dict[ thisMod ] = re.sub(r" ", r"", expr_dict[ thisMod ] )
    expr_dict[ thisMod ] = expr_dict[ thisMod ].replace( '+', 'and' )
    expr_dict[ thisMod ] = expr_dict[ thisMod ].replace( '-', 'or' )
    expr_dict[ thisMod ] = expr_dict[ thisMod ].replace( ',', 'or' )
    expr_dict[ thisMod ] = expr_dict[ thisMod ].strip().replace( '  ', ' ' )
    expr_dict[ thisMod ] = re.sub(r"(\d)(K)", r"\1 \2", expr_dict[ thisMod ] )
    expr_dict[ thisMod ] = re.sub(r"(\d)(\()", r"\1 \2", expr_dict[ thisMod ] )
    expr_dict[ thisMod ] = re.sub(r"\)K", ") K", expr_dict[ thisMod] )
    expr_dict[ thisMod ] = re.sub(r"\)\(", ") (", expr_dict[ thisMod] )
    expr_dict[ thisMod ] = expr_dict[ thisMod ].strip().replace( ' ', 'and' )
    expr_dict[ thisMod ] = expr_dict[ thisMod ].replace( 'and', ' and ' )
    expr_dict[ thisMod ] = expr_dict[ thisMod ].replace( 'or', ' or ' )

# Saved processed pr_expr_dict in module_exprs

mod_to_gene_dict = {}
for thisMod in expr_dict:
    mod_to_gene_dict[thisMod] = re.findall('(K\d{5})', expr_dict[thisMod])

all_genes = list(set(unlistify(list(mod_to_gene_dict.values()))))

first_reactions_dict = {this_gene: [] for this_gene in all_genes}
int_reactions_dict = {this_gene: [] for this_gene in all_genes}
last_reactions_dict = {this_gene: [] for this_gene in all_genes}

from tqdm import tqdm
for this_gene in tqdm(all_genes):
    for thisMod in mod_to_gene_dict:
        if this_gene in mod_to_gene_dict[thisMod][:3]:
            first_reactions_dict[this_gene].append(thisMod)
        elif this_gene in mod_to_gene_dict[thisMod][-3:]:
            last_reactions_dict[this_gene].append(thisMod)
        elif this_gene in mod_to_gene_dict[thisMod]:
            int_reactions_dict[this_gene].append(thisMod)

gene_categories = {}
for this_gene in tqdm(all_genes):
    occs = list(map(len, [first_reactions_dict[this_gene], 
                         int_reactions_dict[this_gene], 
                         last_reactions_dict[this_gene]]))
    gene_categories[this_gene] = occs.index(max(occs))

from ete3 import Tree
import numpy as np
import pandas as pd
from uniqify import uniqify
from unlistify import unlistify
from copy import copy
import pickle
import re
from load_kegg import *

# Getting the set of bacterial abbreviations for each organism in KEGG.
orgNames = []
with open('endo_removed_prok_abbr_kegg.txt', 'r') as f:
    for thisLine in f.readlines():
        orgNames.append( thisLine.strip() )


geneDict = {thisOrg : 
            set(np.genfromtxt('organism_kogenes/' + thisOrg + '.txt').astype(int))
            for thisOrg in orgNames}
gene_ids = sorted( uniqify( unlistify( list( geneDict.values() ) ) ) )

# Getting the tree with internal nodes from gainLoss' ancestral reconstruction.
# Internal nodes are labeled with 'Nx' where x is a number, the root being
# '[N1]' and then the numbers increase.
gl_tree = Tree( 'full_proks_gainLoss_results/TheTree.INodes.ph', format=1 )
ind_tree = Tree('ancs_fastanc_indscores.txt', format = 1)

# # To print the tree.
# print(gl_tree.get_ascii(show_internal=True))

# Getting the root of the tree.
root = gl_tree&'[N1]'
nodes = list( root.traverse() )

gain_positions, loss_positions = [], []
for thisNode in tqdm(nodes):
    for thisChild in thisNode.children:
        gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
        gain_kos = ['K' + str(0) * (5 - len(str(g))) + str(g) for g in gainGenes]
        loss_kos = ['K' + str(0) * (5 - len(str(g))) + str(g) for g in lostGenes]

        gain_positions.append([])
        loss_positions.append([])
        for this_gene in gain_kos:
            if this_gene in gene_categories:
                gain_positions[-1].append(gene_categories[this_gene])

        for this_gene in loss_kos:
            if this_gene in gene_categories:
                loss_positions[-1].append(gene_categories[this_gene])                        

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
            all_fracs.append(0)
    return np.mean(all_fracs)

mean_gain_tallies = {pos: meanAtPos(frac_gain_tallies, pos)
                    for pos in [0, 1, 2]}
mean_loss_tallies = {pos: meanAtPos(frac_loss_tallies, pos)
                    for pos in [0, 1, 2]}

gains, losses = [], []
for thisNode in tqdm(nodes):
    for thisChild in thisNode.children:
        gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
        gain_kos = ['K' + str(0) * (5 - len(str(g))) + str(g) for g in gainGenes]
        loss_kos = ['K' + str(0) * (5 - len(str(g))) + str(g) for g in lostGenes]

        gains +=  gain_kos
        losses += loss_kos

gains, losses = set(gains), set(losses)
ex_gains, ex_losses = gains - losses, losses - gains
ex_gains = set([int(e[1:]) for e in ex_gains])
ex_losses = set([int(e[1:]) for e in ex_losses])

gain_positions, loss_positions = [], []
for this_gene in gains:
    if this_gene in gene_categories:
        gain_positions.append(gene_categories[this_gene])

for this_gene in losses:
    if this_gene in gene_categories:
        loss_positions.append(gene_categories[this_gene])

from collections import Counter

gains_tallies = [Counter(gain_positions[i]) for i in range(len(gain_positions)) 
                 if gain_positions[i]]
loss_tallies = [Counter(loss_positions[i]) for i in range(len(loss_positions))
                if loss_positions[i]]
