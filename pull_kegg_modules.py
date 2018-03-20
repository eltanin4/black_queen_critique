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

# Clearing the expr_dict
import re
expr_dict = pickle.load( open( 'module_exprs/expr_dict.dat', 'rb' ) )
old_expr_dict = deepcopy( expr_dict )
for thisMod in expr_dict:
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
