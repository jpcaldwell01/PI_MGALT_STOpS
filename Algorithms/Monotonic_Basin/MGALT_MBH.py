# FORM: [eval_info,selected] = ...
#       MGALT_MBH(BOD,CONST,OPT,VAR,selected,num_mig,num_isl,count_MBH)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function is the Multiple Gravity-Assist Low-Thrust (MGALT) 
# |     adaptation of the Monotonic Basin Hopping (MBH) function. 
# |
# |     -This function is the main wrapper for the MBH island. Every 
# |     instance of the MBH Island object is used as an input.
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -BOD                (1,1)       [struct]        [unitless]
# |         A struct containing information pertaining to the planetary
# |         bodies. Contains list of bodies, launch windows and ToF, and 
# |         planetary R/V/JD vectors. This struct has dynamic fields and 
# |         will adapt to contain only the necesary information
# |     -CONST              (1,1)       [struct]        [unitless]
# |         A struct containing constants used in the calcs. Contains
# |         values for AU, TU, Sun (rad/mu/rp) and (rad/mu/rp/SOI/per) 
# |         for any bodies used in the optimization scheme. This is a 
# |         dynamic struct and will adapt to contain only the necesary 
# |         information
# |     -OPT                (1,1)       [struct]        [unitless]
# |         A struct containing constants user options. Contains the save 
# |         folder, ToF values, and more structs containing informaiton 
# |         for the island model, cost parameters, weighting parameters, 
# |         and all of the islands used in the optimization process
# |     -VAR                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |     -selected           (Nmig,Nisl) [struct]        [unitless]
# |         Shared solutions from each island for each migration
# |     -num_mig            (1,1)       [int]           [unitless]
# |         The current migration number
# |     -num_isl            (1,1)       [int]           [unitless]
# |         The current island number
# |     -count_MBH          (1,1)       [int]           [unitless]
# |         The current MBH island number
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -eval_info          (1,1)       [struct]        [unitless]
# |         Best solutions, best costs, number of iterations, etc...
# |     -selected           (Nmig,Nisl) [struct]        [unitless]
# |         Shared solutions from each island for each migration
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np
from Other_Functions.genPopn import genPopn
from Island_Model.ISL_modReplacement import ISL_modReplacement
from Algorithms.Monotonic_Basin.MGALT_MBH_function import MGALT_MBH_function
from Island_Model.ISL_bestSoFarCheck import ISL_bestSoFarCheck
from Algorithms.Monotonic_Basin.MGALT_MBH_cluster import MGALT_MBH_cluster
from Island_Model.ISL_modSelection import ISL_modSelection
from IPython import get_ipython
from Other_Functions.MGALT_bestLegs import MGALT_bestLegs
from Other_Functions.randomNum import randomNum


def MGALT_MBH(bod,const,opts,var,selected,num_mig,num_isl,count_MBH):

    ## Pre-Allocation - Loop 1
    

    # Solution Stagnation
    stagnation_first = 0
    
    Nkeep = 100 # Unless you have a buttload of RAM, keep this. (Even if you do, keep it anyways)
    if Nkeep > opts['MBH']['N1_Outer']:
        Nkeep = opts['MBH']['N1_Outer']
    # Make array of 0's for the best solution and cost
    optimal_soln_first = np.zeros((Nkeep,len(var['bin'])))
    
    dates = np.zeros((Nkeep,2*var['transfers']))
    memberDates = np.zeros((1,2*var['transfers']))


    f_best_first = np.full((Nkeep,1),999999999999999999999.)
    
    ## Generate Initial Population - Loop 1
    
    popn_first = genPopn('MBH',opts,opts['MBH'],var)[0]
    
    ## Substitution of best legs
    ps = var['ps']
    if opts['best_Leg_Tracking'] == 'y':
    
        I = np.argsort(popn_first[0,:])[::-1] # Shortcut to get indexes of JD regardless of method
        
             # Transfer 1 -> n-1
        for i1 in reversed(range(1,var['transfers'])):
            for i2 in range(np.shape(popn_first)[0]):
                if np.random.rand(1) <= ps:
                    try: # If MBH is the first algorithm run, this will break because there are no best legs.
                        popn_first[i2,I[i1]:I[i1-1]] = var['bestLegMembers']['leg%i' %(i1+1-var['transfers'])][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(i1+1-var['transfers'])])[0]-1,'int'),I[i1]:I[i1-1]]
                    except:
                        continue
        for i2 in range(np.shape(popn_first)[0]):            
            if np.random.rand(1) <= ps:
                popn_first[i2,I[0]:] = var['bestLegMembers']['leg%i' % (var['transfers']-1)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1)])[0]-1,'int'),I[0]:]
    # popn_first = MGALT_constraintPruning(bod,opts,const,var,popn_first)
    
    ## Add Shared Solutions if Migration Has Occurred
    
    if num_mig > 0:
       popn_first = ISL_modReplacement(opts['island'],selected,num_mig-1,num_isl,popn_first)
    
    ## First Loop
    
    # Progress Display
    get_ipython().magic('clear')
    print('~~~~~~~~~~~ Migrations Completed ',num_mig,'/',opts['island']['Nmig'],'~~~~~~~~~~~')
    print('    Island: ',num_isl+1,'/',opts['island']['Nisl'])
    print('    Algorithm: ',opts['island']['isl_list'][num_isl])
    print('    Loop #1: Searching Global')
    
    # Pass into function
    [f_first,popn_first,is_viable_first,nfeval_first,I,all_fs_trans] = MGALT_MBH_function(bod,const,opts,opts['MBH'],var,popn_first,num_mig,count_MBH,1,Nkeep)
    var = MGALT_bestLegs(var,popn_first,all_fs_trans)

    # Get dates to see clusters for cluster pruning          
    for it in range(Nkeep):
        memberDates[0,:] = popn_first[it,I]       
        dates[it,:] = 1*memberDates
        
    ## Generate Statistics for Current Generation
    avgcost_first = np.mean(f_first)
    mincost_first = min(f_first)
    maxcost_first = max(f_first)
    
    ## If Clusters Were Detected
    
    # Find which rows are viable
    rows = np.nonzero(is_viable_first)[0]

    #If no clusters were found
    if len(rows) == 0:
        
        # Organize Best Solutions So Far
        [f_best_first,optimal_soln_first,stagnation_first] = ISL_bestSoFarCheck(f_first,popn_first,f_best_first,optimal_soln_first,stagnation_first)
        
        # Selecting Solutions For Sharing
        [selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best_first,optimal_soln_first,selected,num_isl,num_mig,opts['island'])
        
        # Eval Info
        eval_info = {}
        eval_info['optimal_soln'] = sorted_optimal_soln 
        eval_info['f_best'] = sorted_f_best
        eval_info['iterations'] = opts['MBH']['N1_Outer']
        eval_info['maxcost'] = maxcost_first 
        eval_info['mincost'] = mincost_first 
        eval_info['avgcost'] = avgcost_first
        eval_info['total_evals'] = nfeval_first
        
        return eval_info,selected,f_first,all_fs_trans
    
    ## Generate Initial Population - Loop 2
    
    popn_second = MGALT_MBH_cluster(bod,opts,opts['MBH'],var,np.hstack((f_first,popn_first)),rows,num_mig)
    
    ## Pre-Allocation - Loop 2
    
    # Solution Stagnation
    stagnation_second = 0
    
    
    # Make array of 0's for the best solution and cost
    optimal_soln_second = np.zeros((np.shape(popn_second)[0],np.shape(popn_second)[1]))
    f_best_second = 999999999999999999999*np.ones((np.shape(popn_second)[0],1))
    
    ## Second Loop
    
    # Progress Display
    print('    Loop #2: Searching Local Basins\n')
    
    # Pass into function
    [f_second,temp,popn_second,nfeval_second] = MGALT_MBH_function(bod,const,opts,opts['MBH'],var,popn_second,num_mig,count_MBH,2)
      
    # Generate Statistics for Current Generation
    avgcost_second = np.mean(f_second)
    mincost_second = min(f_second)
    maxcost_second = max(f_second)
    
    ## Solutions
    
    # Organize Best So Far Solutions
    [f_best_final,optimal_soln_final,stagnation_second] = ISL_bestSoFarCheck(np.vstack(f_first,f_second),np.vstack(popn_first,popn_second),np.vstack(f_best_first,f_best_second),np.vstack(optimal_soln_first,optimal_soln_second),stagnation_second)
    
    # Selecting Solutions For Sharing
    [selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best_final,optimal_soln_final,selected,num_isl,num_mig,opts['island'])
    
    # Eval Info
    eval_info = {}
    eval_info['optimal_soln'] = sorted_optimal_soln
    eval_info['f_best'] = sorted_f_best
    eval_info['iterations'] = opts['MBH']['N1_Outer']+np.shape(popn_second)[0]
    eval_info['maxcost'] = max([maxcost_first,maxcost_second]) 
    eval_info['mincost'] = min([mincost_first,mincost_second]) 
    eval_info['avgcost'] = np.mean([avgcost_first,avgcost_second])
    eval_info['total_evals'] = nfeval_first+nfeval_second
    
    return eval_info,selected,np.vstack((f_first,f_second)),all_fs_trans
