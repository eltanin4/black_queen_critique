import numpy as np

def giveFirstReactions(rxnMat, prodMat, seedVec, sumVec, maxSteps=2):
    """
    Takes in the stoichiometric matrix in the form of a reaction 
    and product matrix, along with a sum vector mentioning the 
    number of reactants in each reaction in the matrix. Also takes
    in a vector of initially seeded metabolites (in seedVec).

    Expands, beginning from the seeded metabolites, a set of reactions
    and metabolites that can be reached iteratively from the former.
    

    RETURNS:
    Returns a set of 'scope-expanded' metabolites and reactions in the 
    stoichiometric matrix form.

    scopeMets is the set of metabolites, with their usual IDs.
    scopeRxns is the set of reactions with my personal IDs 

    NOTE:
    Convert to KEGG IDs before using, please! Else, the apocalypse will
    surely arrive.
    """
    currScopeMets = np.copy(seedVec)
    scopeSize = sum(currScopeMets)
    prevScopeSize = 0.0
    numSteps = 0
    rxnProc = np.zeros(len(prodMat))
    
    while (scopeSize > prevScopeSize and
           numSteps < maxSteps):
        prevScopeSize = scopeSize

        # This generates the rho vector.
        rxnProc = ((np.dot(rxnMat, currScopeMets) - sumVec) == 0) * 1
        # This generates the phi vector.
        currScopeMets = (np.dot(np.transpose(prodMat), rxnProc) + currScopeMets > 0) * 1

        # Recalculating scope size.
        scopeSize = sum(currScopeMets)
        numSteps += 1

    return currScopeMets, rxnProc
