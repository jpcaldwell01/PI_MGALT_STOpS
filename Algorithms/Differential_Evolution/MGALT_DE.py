# FORM: [eval_info,selected] = ...
#       MGALT_DE(BOD,CONST,OPT,VAR,selected,num_mig,num_isl,count_DE)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function is the Multiple Gravity-Assist Low-Thrust (MGALT) 
# |     adaptation of the Differential Evolution (DE) function. 
# |
# |     -This function is the main wrapper for the DE island. Every 
# |     instance of the DE Island object is used as an input.
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
# |     -count_DE           (1,1)       [int]           [unitless]
# |         The current DE island number
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
from Solvers.Indirect import MGALT_IN_FBSM_2D
from Solvers.Direct import MGALT_DIR_FBSM_2D
from Island_Model.ISL_bestSoFarCheck import ISL_bestSoFarCheck
from Algorithms.Differential_Evolution.MGALT_DE_nextGeneration import MGALT_DE_nextGeneration
from Island_Model.ISL_modSelection import ISL_modSelection
from IPython import get_ipython
from joblib import Parallel, delayed
from Other_Functions.clusterPruning import clusterPruning
from Other_Functions.MGALT_bestLegs import MGALT_bestLegs


def MGALT_DE(bod,const,opts,var,selected,num_mig,num_isl,count_DE):

    ## Pre-Allocation
    transfers = len(bod['bodies'])-1
    dates = np.zeros((opts['DE']['Ngen']*opts['DE']['Npop'],2*transfers))
    memberDates = np.zeros((1,2*transfers))
    all_fs = []
    all_fs_trans = np.ndarray((0,var['transfers']))
    
    f = np.zeros((opts['DE']['Npop'],1))
    f_trans = np.zeros((opts['DE']['Npop'],var['transfers'])) 
    avgcost = np.zeros((opts['DE']['Ngen'],1))
    mincost = np.zeros((opts['DE']['Ngen'],1))
    maxcost = np.zeros((opts['DE']['Ngen'],1))
    
    # Solver Selection
    
    if opts['solver'] == 'MGALT_IN_FBSM_2D':
        feval = getattr(MGALT_IN_FBSM_2D, 'MGALT_IN_FBSM_2D')
        I = [0] + [11*x-1 for x in range(1,var['transfers'])] + [-1] + [11*x for x in range(1,var['transfers'])] # Cluster pruning indexes

    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
        feval = getattr(MGALT_DIR_FBSM_2D, 'MGALT_DIR_FBSM_2D')
        I = [0] + [(opts['thrust']['Nseg']*2+5)*x-1 for x in range(1,var['transfers'])] + [-1] + [(opts['thrust']['Nseg']*2+5)*x for x in range(1,var['transfers'])] # Cluster pruning indexes

    else:
        raise Exception('Invalid Solver Selection')
    # Solution Stagnation
    stagnation = 0 
    
    # Number of evals
    nfeval = 0
    
    # Make arrays of 0's for the best solution and cost
    optimal_soln = np.zeros((opts['DE']['Npop'],len(var['bin'])))
    f_best = 999999999999999999999*np.ones((opts['DE']['Npop'],1))
    
    
    
    ## Generate Initial Population
    
    popn = genPopn('DE',opts,opts['DE'],var)[0]
    
    ## Add Shared Solutions if Migration Has Occurred
    
    if num_mig > 0:
        popn = ISL_modReplacement(opts['island'],selected,num_mig-1,num_isl,popn)
    
    ## Iterate Generations
    
    # Progress Display
    get_ipython().magic('clear')
    print('~~~~~~~~~~~ Migrations Completed ',num_mig,'/',opts['island']['Nmig'],'~~~~~~~~~~~')
    print('    Island: ',num_isl+1,'/',opts['island']['Nisl'])
    print('    Algorithm: ',opts['island']['isl_list'][num_isl]) 
    
    # Calculate Fitness Value for Each Member
    if opts['parallel'] == 'y':
        
        # Preallocate
        # data(1:OPT.DE(count_DE).Npop) = parallel.FevalFuture
        # fh = str2func(string(OPT.solver))      # https://www.mathworks.com/help/matlab/ref/str2func.html
        
        # Display Info
        ll_par = print('\n    Parallel Processing   \n')
        ll_RAM = print('    Watch CPU/RAM usage   \n\n')
        ll_mem = print('\r','    Generating ',opts['DE']['Npop'],' Members',end='')
        
        # Cost
        f,plot_vars,f_trans = zip(*Parallel(n_jobs=8)(delayed(feval)(popn[member],bod,const,opts,var) for member in range(opts['DE']['Npop'])))
        f = np.asarray(f)
        if opts['best_Leg_Tracking'] == 'y':
            f_trans = np.squeeze(np.asarray(f_trans))
            var = MGALT_bestLegs(var,popn,np.squeeze(f_trans))

         
        nfeval = nfeval + opts['DE']['Npop']
    
        # Collect results as they become available
        
    elif opts['parallel'] == 'n':
    
        # Display Info
        # ll_mem = fprintf('    Member: 1/#1.0f\n',OPT.DE(count_DE).Npop)
        
        # Cost
        for member in range(opts['DE']['Npop']):
            
            # Generation Display
            # fprintf(repmat('\b',1,ll_mem))
            ll_mem = print('\r', '   Member: ',str(member+1),'/',str(opts['DE']['Npop']), end = '')  
    
            # Cost 
            f[member],plot_vars,f_trans[member,:] = feval(popn[member],bod,const,opts,var)
            
        if opts['best_Leg_Tracking'] == 'y':
            var = MGALT_bestLegs(var,popn,f_trans)
            
        nfeval = nfeval + opts['DE']['Npop']

    else:
        raise Exception("Invalid Parallel Processing Input. Must be 'y' or 'n'")
    
    # Progress Display
    # print('~~~~~~~~~~~ Migrations Completed ',num_mig,'/',opts['island']['Nmig'],'~~~~~~~~~~~')
    # print('    Island: ',num_isl+1,'/',opts['island']['Nisl'])
    # print('    Algorithm: ',opts['island']['isl_list'][num_isl]) 
    
    for gen in range(opts['DE']['Ngen']):
        
        # Generate Statistics for Current Generation
        avgcost[gen] = np.mean(f) 
        mincost[gen] = min(f) 
        maxcost[gen] = max(f)
       
        # Organize Best So Far Solutions
        [f_best,optimal_soln,stagnation] = ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation)
       
        # Progress Display
        # fprintf(repmat('\b',1,ll_gen))
        ll_gen = print("\r",'   Generation: ',str(gen+1),'/',str(opts['DE']['Ngen']),end='       ')
       
        # Get dates to see clusters            
        for it in range(np.shape(popn)[0]):
            memberDates[0,:] = popn[it,I]       
            dates[(gen)*opts['DE']['Npop'] + it,:] = 1*memberDates
            
        # Select Mating Pool Mate to Create Next Generation
        [popn,f,nfeval,f_trans] = MGALT_DE_nextGeneration(bod,const,opts,opts['DE'],var,optimal_soln,f_best,nfeval,count_DE)
        if opts['best_Leg_Tracking'] == 'y':
            var = MGALT_bestLegs(var,popn,np.squeeze(f_trans))
            all_fs_trans = np.append(all_fs_trans,np.squeeze(f_trans),axis=0)


        all_fs = np.append(all_fs,f)
    
    ## Selecting Solutions For Sharing
    [f_best,optimal_soln,stagnation] = ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation)
    [selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best,optimal_soln,selected,num_isl,num_mig,opts['island'])
    
    
    
    ## Eval Info

    eval_info = {}
    eval_info['optimal_soln'] = sorted_optimal_soln
    eval_info['f_best'] = sorted_f_best
    eval_info['iterations'] = opts['DE']['Ngen']
    eval_info['maxcost'] = maxcost 
    eval_info['mincost'] = mincost 
    eval_info['avgcost'] = avgcost
    eval_info['total_evals'] = nfeval
    
    plot_me = 0 # Change to 1 to show cluster
    var = clusterPruning(dates,opts,var,'DE',num_mig,I,plot_me)
    
    return eval_info,selected,all_fs,all_fs_trans
