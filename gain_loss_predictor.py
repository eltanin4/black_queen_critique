from give_scope import giveScope

def giveSeedVec(sources):
    # Generating a vector of all metabolites initially provided, i.e. seeds.
    seedVec = np.zeros(len(rxnMat.T))
    # media = list(sources['kegg_id'].values)
    # media = list(sources)
    media = list(sources & set(cpds))
    seedVec[[kegg_to_id[e] for e in Currency + media]] = 1
    return seedVec


# Instantiating default rich medium.
seeds_df = pd.read_csv('seeds_from_vitkup.csv')
default_media = set(list(seeds_df['kegg_id'].values))
defaultSeedVec = giveSeedVec(default_media)

#--------------------------------------------------------------------------------------------------

def give_number( string ):
    return int( string[-5:] )

#--------------------------------------------------------------------------------------------------

def genotypeToRxns(genotype):
    kos = [pos_to_gene_map[i] for i in np.where((genotype > 0.5) * 1)[0]]
    refDF = pd.read_csv('KOREF.txt', delimiter='\t', sep='delimiter', header=None, names=['id', 'rid'])
    refDF['id'] = refDF['id'].apply( give_number )
    refDF['rid'] = refDF['rid'].apply( give_number )
    can = refDF.loc[refDF['id'].isin(kos),:]['rid'].values
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in can]
    can = list((set(can)) & set(rxns))
    return can

#--------------------------------------------------------------------------------------------------

def inferProducible(genotype, sources=defaultSeedVec):
    """
    Given a genotype, infers what reactions it can perform, and then calculates 
    all producible metabolites in a given environment. 

    If no environment is supplied, assumes a rich environment with all single 
    carbon source enrichments.
    """
    # Getting all reactions that the genotype can perform.
    can = genotypeToRxns(genotype)

    # Generating the seed vector from the environmental sources.
    currSeedVec = giveSeedVec(sources)

    avrxns = [rxn_kegg_to_id[e] for e in can]
    scopeMets = giveScope(rxnMat[avrxns], prodMat[avrxns], currSeedVec, sumRxnVec[avrxns])[0]
    return set([id_to_kegg[e] for e in np.where(scopeMets)[0]])

#--------------------------------------------------------------------------------------------------

def coreFromReactions(sources, reactionSet):
    avrxns = [rxn_kegg_to_id[e] for e in reactionSet]
    scopeMets, avScopeRxn = giveScope(rxnMat[avrxns], prodMat[avrxns], giveSeedVec(sources), 
                                     sumRxnVec[avrxns])
    scopeRxns = np.nonzero(avrxns * avScopeRxn)[0]
    return set([id_to_kegg[e] for e in set(np.where(scopeMets)[0])]) & set(Core)

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

# Current environment is E, partner genotype is P.
# Ancestral (focal) microbe is anc.
ancNode = gl_tree&(tip_to_name_segata_dict[kegg_to_tip_name_dict['elc']])
anc = reconAncestor(anc_recon_table, ancNode)
lossDepsDict = {}

for prtID in tqdm(range(len(orgNames))):
    for envID in tqdm(range(len(default_media))):
        # Specify partner here.
        thisOrg = orgNames[prtID]
        prtNode = gl_tree&(tip_to_name_segata_dict[kegg_to_tip_name_dict[thisOrg]])
        P = reconAncestor(anc_recon_table, prtNode)

        # Specify environment here.
        E = set([list(default_media)[envID]])

        # Calculating what P can secrete/produce in E. Currently interchangeable.
        beta = inferProducible(P, E)

        # Now calculating producible cores.
        # LEGEND:   n1 is anc on E alone.
        #           n2 is anc on E + P.
        #           n3 is anc on P alone.
        #           n4 is P alone.
        n1 = giveCoresProduced(E, anc)
        n2 = giveCoresProduced(E | beta, anc)
        n3 = giveCoresProduced(beta, anc)
        n4 = giveCoresProduced(E, P)

        # Calculating candidates for dependency evolution.
        # Checking for cases of possible gene loss first, 
        # also, can serve as a sanity check.
        lossCores = n3 & n1

        # Getting the set of all KEGG pathway modules.
        pr_expr_dict = pickle.load( open( 'module_exprs/pr_expr_dict.dat', 'rb' ) )
        m2rxn = pd.read_csv('module_to_reaction.txt', delimiter='\t', sep='delimiter', header=None, names=['mid', 'rid'])

        # First getting the reactions that anc can perform.
        rAnc = genotypeToRxns(anc)

        for thisMod in tqdm(pr_expr_dict):
            # Initializing the module dictionary.
            lossDepsDict[thisMod] = []

            # Getting all genes in this module.
            rxnsInMod = m2rxn.loc[m2rxn['mid'] == 'md:' + thisMod,:]['rid'].values
            rxnsInMod = [m[3:] for m in rxnsInMod]
            modPresent = list((set(rxnsInMod)) & set(rxns))
            modPresent = [rxn_kegg_to_id[e] for e in modPresent]

            # Creating a mutant knockout genotype of anc, say mAnc, in terms of reactions only.
            mAnc = [e for e in rAnc if e not in modPresent]

            # Getting the produciblity post loss.
            nmut = coreFromReactions(E, mAnc)

            # Checking which candidate cores cannot be produced in mAnc.
            nDeps = [thisCore for thisCore in lossCores 
                     if thisCore not in nmut]

            # Updating a global dictionary.
            lossDepsDict[thisMod] += nDeps

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

# Using phylogenetic history to make guesses about gains and losses involved in
# dependency evolution.
lossDepsList, cglDepsList = [], []
tlossDepsList, tcglDepsList = [], []
E = unuts.copy()
for thisNode in tqdm(nodes):
    for thisChild in thisNode.children:
        # Defining the ancestor and descendant.
        anc = thisNode.genotype
        rAnc = genotypeToRxns(anc)

        # Figuring out what the ancestor can do.
        n1 = giveCoresProduced(E, anc)
        n2 = giveCoresProduced(E | ubyps, anc)
        n3 = giveCoresProduced(ubyps, anc)

        # Cataloging the gains and losses along the branch.
        gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)

        # First checking for just-loss dependencies.
        mAnc = np.copy(anc)
        lostGenePos = [gene_to_pos_map[i] for i in lostGenes]
        mAnc[lostGenePos] = 0
        n4 = giveCoresProduced(E, mAnc)
        n5 = giveCoresProduced(ubyps, mAnc)

        lossDeps = (n1 & n5) & set(Core).difference(n4)
        lossDepsList.append(lossDeps)

        # Now checking for CGL (conjugate gains and losses) mediated dependencies.
        gAnc = np.copy(anc)
        gainGenePos = [gene_to_pos_map[i] for i in gainGenes]
        gAnc[gainGenePos] = 1
        n6 = giveCoresProduced(E, gAnc)
        n7 = giveCoresProduced(ubyps, gAnc)
        n8 = giveCoresProduced(E, thisChild.genotype)
        n9 = giveCoresProduced(ubyps, thisChild.genotype)

        # Currently CGLs operate under the stringent criterion, i.e. gains allow new core production.
        cglDeps = (n1 & (n7 & set(Core).difference(n5))) & set(Core).difference(n8)
        cglDepsList.append(cglDeps)

        # Sifting out terminal (extant) genotypes.
        if not thisChild.children:
            tlossDepsList.append(lossDeps)
            tcglDepsList.append(cglDeps)

# Now sorting, filtering and analyzing dependency differences.
numLosses = np.array([len(e) for e in tlossDepsList])
numCGLs = np.array([len(e) for e in tcglDepsList])
eventRatios = np.divide(numLosses, numCGLs)

plt.hist(eventRatios, bins=10)
plt.show()

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

# What about phylogenetic histories where a gain on a previous branch 
# is correlated with a loss in the subsequent branch? Those are also 
# technically CGL events, and this snippet searches for the events consistent
# with this. (Note: only searches for adjacent branches, ignores distantly 
# related events).
pcglDepsList = []
ptcglDepsList = []
E = unuts.copy()
for thisNode in tqdm(nodes):
    for thisChild in thisNode.children:
        if thisNode.up:
            # Defining the ancestor and descendant.
            anc = thisNode.genotype
            rAnc = genotypeToRxns(anc)
            parAnc =  thisNode.up.genotype


            # Figuring out what the parent of the ancestor can do.
            s1 = giveCoresProduced(E, parAnc)
            s2 = giveCoresProduced(ubyps, parAnc)

            # Cataloging the gains and losses along the branch.
            gainGenes, lostGenes = giveGainsAndLosses(thisNode, thisChild)
            pgainGenes, plostGenes = giveGainsAndLosses(thisNode.up, thisNode)

            # Now checking for CGL (conjugate gains and losses) mediated dependencies.
            gAnc = np.copy(parAnc)
            gainGenePos = [gene_to_pos_map[i] for i in pgainGenes]
            gAnc[gainGenePos] = 1
            s3 = giveCoresProduced(ubyps, gAnc)
            
            # Now performing the subsequent gene losses.
            glAnc = np.copy(gAnc)
            lossGenePos = [gene_to_pos_map[i] for i in lostGenes]
            glAnc[lossGenePos] = 0
            s4 = giveCoresProduced(E, glAnc)
            
            pcglDeps = (s1 & (s3 & set(Core).difference(s2))) & set(Core).difference(s4)
            pcglDepsList.append(pcglDeps)

            # Sifting out terminal (extant) genotypes.
            if not thisChild.children:
                ptcglDepsList.append(pcglDeps)
