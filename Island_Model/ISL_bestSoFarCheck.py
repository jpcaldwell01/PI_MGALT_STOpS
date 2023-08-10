# FORM: [f_best,optimal_soln,stagnation] = ...
#       ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation)
#
# |-----------------------------------------------------------------------
# | NOTES:
# |     -Best So Far Check funciton. Sorts the solutions into the best
# |     order
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -f                  (Npop,1)  	[float]      	[unitless]
# |         The cost associated with every variable string in the 
# |         population
# |     -popn               (Npop,Nvar) [float]         [unitless]
# |         The current population of Npop variable strings
# |     -f_best             (Npop,1)	[float]         [unitless]
# |         The cost associated with every variable string in the 
# |         population in best to worst order
# |     -optimal_soln    	(Npop,Nvar) [float]         [unitless]
# |         New array of optimal variables strings associated with costs 
# |         in 'f_best'
# |     -stagnation       	(1,1)       [int]           [unitless]
# |         Measure of if the algorithm is stagnating. Counts how many 
# |         generations in a row with no better solution
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -f_best             (Npop,1)  	[float]      	[unitless]
# |         New array of best costs in order
# |     -optimal_soln    	(Npop,Nvar) [float]         [unitless]
# |         New array of optimal variables strings associated with costs 
# |         in 'f_best'
# |     -stagnation       	(1,1)       [int]           [unitless]
# |         Measure of if the algorithm is stagnating. Counts how many 
# |         generations in a row with no better solution
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np

def ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation):

    ## Setup
    
    Npop = len(f)
    found_better = 0
    
    for j in range(Npop):
        
        if len(f) == 0:
            break
       
        while any(min(f) == f_best):
            cut = np.argmin(f)
            f = np.delete(f,cut) 
            popn = np.delete(popn,cut,0)
          
            if len(f) == 0:
                break
       
        if len(f) == 0:
            break
       
        if any(min(f) < f_best): # If it's better than any of them
            found_better = 1
            replace_loc = np.argmax(f_best)
            f_best[replace_loc] = min(f)
            ind = np.argmin(f)
            optimal_soln[replace_loc] = popn[ind]
            f = np.delete(f,ind) 
            popn = np.delete(popn,ind,0)
       
    
    if found_better:
        stagnation = 0
    else:
        stagnation = stagnation + 1 # If a better tour was not found
    
    return f_best,optimal_soln,stagnation