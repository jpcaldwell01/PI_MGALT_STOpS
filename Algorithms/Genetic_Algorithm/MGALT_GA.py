# FORM: [eval_info,selected] = ...
#       MGALT_GA(BOD,CONST,OPT,VAR,selected,num_mig,num_isl,count_GA)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function is the Multiple Gravity-Assist Low-Thrust (MGALT) 
# |     adaptation of the Genetic Algorithm (GA) function. 
# |
# |     -This function is the main wrapper for the GA island. Every 
# |     instance of the GA Island object is used as an input.
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
# |     -count_GA           (1,1)       [int]           [unitless]
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
from Island_Model.ISL_modSelection import ISL_modSelection
from Solvers.Indirect import MGALT_IN_FBSM_2D
from Solvers.Direct import MGALT_DIR_FBSM_2D
from Island_Model.ISL_bestSoFarCheck import ISL_bestSoFarCheck
from Algorithms.Genetic_Algorithm.MGALT_GA_nextGeneration import MGALT_GA_nextGeneration
from IPython import get_ipython
from joblib import Parallel, delayed
from Other_Functions.MGALT_bestLegs import MGALT_bestLegs

def MGALT_GA(bod,const,opts,var,selected,num_mig,num_isl,count_GA):

    #%% Pre-Allocation
    transfers = len(bod['bodies'])-1
    all_fs = []
    all_fs_trans = np.ndarray((0,var['transfers']))
    dates = np.zeros((opts['GA']['Ngen']*opts['GA']['Npop'],2*transfers))
    memberDates = np.zeros((1,2*transfers))
    # Solver Selection
    
    if opts['solver'] == 'MGALT_IN_FBSM_2D':
        feval = getattr(MGALT_IN_FBSM_2D, 'MGALT_IN_FBSM_2D')
        I = [0] + [11*x-1 for x in range(1,var['transfers'])] + [-1] + [11*x for x in range(1,var['transfers'])] # Cluster pruning indexes

    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
        feval = getattr(MGALT_DIR_FBSM_2D, 'MGALT_DIR_FBSM_2D')
        I = [0] + [(opts['thrust']['Nseg']*2+5)*x-1 for x in range(1,var['transfers'])] + [-1] + [(opts['thrust']['Nseg']*2+5)*x for x in range(1,var['transfers'])] # Cluster pruning indexes

    else:
        raise Exception('Invalid Solver Selection')
    
    # Cost Parameters
            ########### MATLAB CODE (note difference in count_GA)###########
            # f = zeros(OPT.GA(count_GA).Npop,1)
            # avgcost = zeros(OPT.GA(count_GA).Ngen,1)
            # mincost = zeros(OPT.GA(count_GA).Ngen,1)
            # maxcost = zeros(OPT.GA(count_GA).Ngen,1)
    f = np.zeros((opts['GA']['Npop'],1)) 
    f_trans = np.zeros((opts['GA']['Npop'],var['transfers']))
    avgcost = np.zeros((opts['GA']['Ngen'],1))
    mincost = np.zeros((opts['GA']['Ngen'],1))
    maxcost = np.zeros((opts['GA']['Ngen'],1))
    
    # Solution Stagnation
    stagnation = 0 
    
    # Number of evals
    nfeval = 0
    # Make arrays of 0's for the best solution and cost
                            ########### MATLAB CODE (note difference in count_GA)###########
                        # optimal_soln = zeros(OPT.GA(count_GA).Npop,size(VAR.bin,2))
                        # f_best = 999999999999999999999*ones(OPT.GA(count_GA).Npop,1)
    optimal_soln = np.zeros((opts['GA']['Npop'],len(var['bin'])))
    f_best =  999999999999999999999*np.ones((opts['GA']['Npop'],1))
    
    
    ## Generate Initial Population
                        # popn = genPopn('GA',OPT,OPT.GA(count_GA),VAR)[0]
    popn = genPopn('GA',opts,opts['GA'],var)[0]
    
    # popn = MGALT_constraintPruning(bod,opts,const,var,popn)
    # all_popn = 1*popn
    
    ## Add Shared Solutions if Migration Has Occurred
    
    if num_mig > 0:
        popn = ISL_modReplacement(opts['island'],selected,num_mig-1,num_isl,popn)
        # all_popn = np.vstack((all_popn,popn))

    
    # Fix Initial Thruster Angles
    
    
    
    #%% Iterate Generations
    
    # Progress Display

    get_ipython().magic('clear')
    print('~~~~~~~~~~~ Migrations Completed ',num_mig,'/',opts['island']['Nmig'],'~~~~~~~~~~~')
    print('    Island: ',num_isl+1,'/',opts['island']['Nisl'])
    print('    Algorithm: ',opts['island']['isl_list'][num_isl])
    
    if opts['parallel'] == 'y':
        ll_par = print('\n    Parallel Processing   \n')
        ll_RAM = print('    Watch CPU/RAM usage   \n\n')
    
    for num_gen in range(opts['GA']['Ngen']):
        
        # Calculate Fitness Value for Each Member
        
        if opts['parallel'] == 'y':
            
            # Display Info
            print('\r','Generation: ',num_gen+1,'/',opts['GA']['Ngen'],end=" ")
            
            # Cost
            f,plot_vars,f_trans = zip(*Parallel(n_jobs=8)(delayed(feval)(popn[member],bod,const,opts,var) for member in range(opts['GA']['Npop'])))
            f = np.asarray(f)
            all_fs = np.append(all_fs,f)
            if opts['best_Leg_Tracking'] == 'y':
                all_fs_trans = np.append(all_fs_trans,np.squeeze(f_trans),axis=0)
                var = MGALT_bestLegs(var,popn,np.squeeze(f_trans))

            
        elif opts['parallel'] == 'n':
        
            # Display Info
            print("\r",'   Generation: ',str(num_gen+1),'/',str(opts['GA']['Ngen']),end='')
                
            for member in range(opts['GA']['Npop']):
                
                # Cost
                f[member],plot_vars,f_trans[member,:] = feval(popn[member],bod,const,opts,var)
            
            all_fs = np.append(all_fs,f)
            if opts['best_Leg_Tracking'] == 'y':
                all_fs_trans = np.append(all_fs_trans,f_trans,axis=0)
                var = MGALT_bestLegs(var,popn,f_trans)
                
        else:
            raise Exception("Invalid Parallel Processing Input. Must be 'y' or 'n'")
            
        nfeval = nfeval + opts['GA']['Npop']
        
        # Generate Statistics for Current Generation
        avgcost[num_gen] = np.mean(f) 
        mincost[num_gen] = min(f) 
        maxcost[num_gen] = max(f)
        
        # Organize Best So Far Solutions
        [f_best,optimal_soln,stagnation] = ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation)
            
        # Get dates to see clusters            
        for it in range(np.shape(popn)[0]):
            memberDates[0,:] = popn[it,I]       
            dates[(num_gen)*opts['GA']['Npop'] + it,:] = 1*memberDates
        
        # Select Mating Pool Mate to Create Next Generation
        popn = MGALT_GA_nextGeneration(bod,opts,opts['GA'],var,popn,f)
        # popn = MGALT_constraintPruning(bod,opts,const,var,popn)
        # all_popn = np.vstack((all_popn,popn))
    
    #%% Selecting Solutions For Sharing
    
    [selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best,optimal_soln,selected,num_isl,num_mig,opts['island'])
    
    ## Eval Info
    eval_info = {}
    eval_info['optimal_soln'] = sorted_optimal_soln
    eval_info['f_best'] = sorted_f_best
    eval_info['iterations'] = opts['GA']['Ngen']
    eval_info['maxcost'] = maxcost 
    eval_info['mincost'] = mincost 
    eval_info['avgcost'] = avgcost
    eval_info['total_evals'] = nfeval
    
    return eval_info,selected,all_fs,all_fs_trans
