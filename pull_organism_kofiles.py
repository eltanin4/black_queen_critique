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

#-------------------------------------------------------------------------
# This code extracts all KEGG KO pages for a given set of organism abbrs.
#-------------------------------------------------------------------------
organism_url_holder = 'http://rest.kegg.jp/link/ko/'

handler = urllib.request.ProxyHandler({'http': 'proxy.ncbs.res.in:3128'})
opener = urllib.request.build_opener(handler)
urllib.request.install_opener(opener)
urllib.request.ProxyHandler({'ftp': 'proxy.ncbs.res.in:3128'})

# Getting the overall save directory.
saveDir = 'organism_konums/'
if not os.path.exists( saveDir ):
    os.mkdir( saveDir )

# Getting the names of the abbrevations from KEGG.
orgNames = []
with open('endo_removed_prok_abbr_kegg.txt', 'r') as f:
    for thisLine in f.readlines():
        orgNames.append( thisLine.strip() )

# Pulling all organism KO numbers.
for orgIndex, orgAbbr in enumerate( orgNames ):
    print_progress_bar( orgIndex, len( orgNames ), 'Pulling all organism KO numbers')
    current_url = organism_url_holder + str( orgAbbr )
    path = saveDir + str( orgAbbr ) + '.html'
    if not os.path.isfile(path):
        try:
            urllib.request.urlretrieve( current_url, path )
        except:
            pass

def find_between_r(s, first, last):
    try:
        start = s.rindex(first) + len(first)
        end = s.rindex(last, start)
        return s[start:end]
    except ValueError:
        return ""

# Converting organism EC number files to pure lists of EC numbers
# to feed to NetCmpt.
for orgIndex, orgAbbr in enumerate( orgNames ):
    print_progress_bar( orgIndex, len( orgNames ), 'Purging all irrelevant content from downloaded files')
    # Get the list of strings, each row corresponding to an EC number.
    try:
        with open( saveDir + str( orgAbbr ) + '.html', 'r' ) as orgFile:
            orgString = orgFile.readlines()
    except:
        continue

    # Extracting the EC number from each line.
    thisOrgECList = []
    for thisLine in orgString:
        thisOrgECList.append( find_between_r( thisLine, 'ko:K', '\n' ) ) 

    with open( saveDir + str( orgAbbr ) + '.txt', 'w' ) as outFile:
        outFile.write( '\n'.join( thisOrgECList ) )

    # Creating file for NetCmpt.
    with open( 'all_prok_konums.txt', 'a' ) as outFile:
        outFile.write( str( orgAbbr ) + ' ' + ' '.join( thisOrgECList ) + '\n' )


# Now getting all reactions in each organism.
# Getting the overall save directory.
rxn_saveDir = 'organism_reactions/'
if not os.path.exists( rxn_saveDir ):
    os.mkdir( rxn_saveDir )

ko_saveDir = 'organism_kogenes/'
if not os.path.exists( ko_saveDir ):
    os.mkdir( ko_saveDir )

# Extracting and saving all KEGG reaction IDs.
num_failed = 0
for orgIndex, orgAbbr in enumerate( orgNames ):
    print_progress_bar( orgIndex, len( orgNames ), 'Saving reaction IDs for each organism')
    try:
        orgECs = pd.read_table( saveDir + str( orgAbbr ) + '.txt', delimiter=' ', sep='delimiter', header=None, names=['id'] )
    except:
        num_failed += 1
    refDF = pd.read_csv('KOREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])

    def give_number( string ):
        return int( string[-5:] )

    refDF['id'] = refDF['id'].apply( give_number )
    refDF['rid'] = refDF['rid'].apply( give_number )
    orgECs.drop_duplicates()
    koGenes = pd.merge( orgECs, refDF, how='left', on='id' ).dropna( )[ 'id' ].tolist( )
    keggRxns = pd.merge( orgECs, refDF, how='left', on='id' ).dropna( )[ 'rid' ].tolist( )
    keggRxns = np.array( [ int( thisRxn ) for thisRxn in keggRxns ] )

    np.savetxt( ko_saveDir + str( orgAbbr ) + '.txt', koGenes )
    np.savetxt( rxn_saveDir + str( orgAbbr ) + '.txt', keggRxns )
