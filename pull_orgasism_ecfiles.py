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
# This code extracts all KEGG webpages in the pathway map 01120.
#-------------------------------------------------------------------------
organism_url_holder = 'http://rest.kegg.jp/link/enzyme/'

handler = urllib.request.ProxyHandler({'http': 'proxy.ncbs.res.in:3128'})
opener = urllib.request.build_opener(handler)
urllib.request.install_opener(opener)
urllib.request.ProxyHandler({'ftp': 'proxy.ncbs.res.in:3128'})

# Getting the overall save directory.
saveDir = 'organism_ecnums/'
if not os.path.exists( saveDir ):
    os.mkdir( saveDir )

# Getting a list of all KEGG prokaryotic abbreviations from file.
with open( 'prok_kegg_abbreviations.txt', 'r' ) as orgFile:
    orgNames = orgFile.readlines()
orgNames = [ x.strip() for x in orgNames ]  

# Pulling all organism EC numbers.
for orgIndex, orgAbbr in enumerate( orgNames ):
    print_progress_bar( orgIndex, len( orgNames ), 'Pulling all organism EC numbers')
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
        thisOrgECList.append( find_between_r( thisLine, 'ec:', '\n' ) ) 

    with open( saveDir + str( orgAbbr ) + '.txt', 'w' ) as outFile:
        outFile.write( '\n'.join( thisOrgECList ) )

    with open( saveDir + str( orgAbbr ) + '.txt', 'w' ) as outFile:
        outFile.write( str( orgAbbr ) + ' ' + ' '.join( thisOrgECList ) )

# Creating file for NetCmpt.
for orgIndex, orgAbbr in enumerate( orgNames ):
    print_progress_bar( orgIndex, len( orgNames ), 'Compiling file for NetCmpt')
    with open( 'all_prok_ecnums.txt', 'a' ) as outFile:
        outFile.write( str( orgAbbr ) + ' ' + ' '.join( thisOrgECList ) + '\n' )


# Now getting all reactions in each organism.
# Getting the overall save directory.
rxn_saveDir = 'organism_reactions/'
if not os.path.exists( rxn_saveDir ):
    os.mkdir( rxn_saveDir )

# Extracting and saving all KEGG reaction IDs.
for orgIndex, orgAbbr in enumerate( orgNames ):
    print_progress_bar( orgIndex, len( orgNames ), 'Saving reaction IDs for each organism')
    try:
        orgECs = pd.read_csv( saveDir + str( orgAbbr ) + '.txt', sep='delimiter', header=None, names=['id'] )
    except:
        continue
    refDF = pd.read_csv('ECREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])

    def numeric_id(string):
        return '.'.join(re.findall('\d+', string))

    orgECs['id'] = orgECs['id'].apply(numeric_id)
    refDF['id'] = refDF['id'].apply(numeric_id)
    orgECs.drop_duplicates()
    orgECs.describe()
    refDF.describe()
    keggRxns = pd.merge( orgECs, refDF, how='left', on='id' ).dropna( )[ 'rid' ].tolist( )
    keggRxns = np.array( [ int( thisRxn[ 4: ] ) for thisRxn in keggRxns ] )

    np.savetxt( rxn_saveDir + str( orgAbbr ) + '.txt', keggRxns )
