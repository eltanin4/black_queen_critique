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
# First getting a hierarchical list of all KEGG organism abbreviations
# and names.
#-------------------------------------------------------------------------
# kegg_dict = {}
# with open( 'kegg_abbr_htext.keg', 'r' ) as modFile:
#     # Going over the htext line-by-line.
#     for thisLine in modFile.readlines():
#         # Checking if new category reached.
#         if thisLine[0] == 'E':
#             abbr = thisLine.split()[ 1 ]
#             kegg_dict[ abbr ] = ' '.join( thisLine.split()[ 2: ] )

# # Saved names as kegg_abbr_to_name_dict.dat
spdf = pd.read_excel( 'gaurav_species_phylo.xlsx' )
# species = [ e.replace( 'str.', '' ) for e in spdf[ 'ncbi_name' ].values ]

kegg_dict = pickle.load( open( 'kegg_abbr_to_name_dict.dat', 'rb' ) ) 
inv_kegg_dict = { value: key for key, value in kegg_dict.items() }

#-------------------------------------------------------------------------
# This code extracts all KEGG webpages for KEGG organisms.
#-------------------------------------------------------------------------
# module_url_holder = 'http://www.kegg.jp/kegg-bin/show_organism?org='

# handler = urllib.request.ProxyHandler({'http': 'proxy.ncbs.res.in:3128'})
# opener = urllib.request.build_opener(handler)
# urllib.request.install_opener(opener)
# urllib.request.ProxyHandler({'ftp': 'proxy.ncbs.res.in:3128'})

# Getting the overall save directory.
saveDir = 'kegg_org_pages/'
if not os.path.exists( saveDir ):
    os.mkdir( saveDir )

# # Pulling all module pages in KEGG MODULE.
# for ti, ts in enumerate( kegg_dict ):
#     print_progress_bar( ti, len( kegg_dict ), 'Getting all organism')
#     current_url = module_url_holder + str( ts )
#     path = saveDir + str( ts ) + '.html'
#     if not os.path.isfile(path):
#         urllib.request.urlretrieve( current_url, path )

#-------------------------------------------------------------------------
# This code gets all NCBI accession numbers from the webpages for each 
# organism downloaded.
#-------------------------------------------------------------------------
# Getting the logical expression string corresponding to each module.
# ncbi_acc_dict = {}
# for ti, ts in enumerate( kegg_dict ):
#     print_progress_bar( ti, len( kegg_dict ), 'Finding all accession numbers')
#     # Get the list of strings, each row corresponding to an EC number.
#     try:
#         with open( saveDir + str( ts ) + '.html', 'r' ) as modFile:
#             modString = BeautifulSoup(modFile, "lxml").text
#     except:
#         continue

#     # Extracting the right expression string.
#     ncbi_acc_dict[ ts ] = find_between_r( modString, 'Assembly: ', ')' )[:15]

# lem and bbev are missing
# Saved the inverse ncbi_to_kegg_dict as ncbi_to_kegg_dict.dat
ncbi_to_kegg_dict = pickle.load( open( 'ncbi_to_kegg_dict.dat', 'rb' ) )
ncbi_to_kegg_df = pd.DataFrame( list( ncbi_to_kegg_dict.items() ), 
                                columns=[ 'ncbi_acc', 'kegg_abbr' ] )
spdf = spdf.merge( ncbi_to_kegg_df, on='ncbi_acc' )
