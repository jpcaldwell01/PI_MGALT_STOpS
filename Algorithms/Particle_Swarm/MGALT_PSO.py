# FORM: [eval_info,selected] = ...
#       MGALT_PSO(solver,inputs,opt_PSO,selected,mig,isl)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function is the Multiple Gravity-Assist Low-Thrust (MGALT) 
# |     adaptation of the Particle Swarm Optimization (PSO) function. 
# |
# |     -This function is the main wrapper for the PSO island. Every 
# |     instance of the PSO Island object is used as an input.
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
# |     -count_PSO          (1,1)       [int]           [unitless]
# |         The current PSO island number
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
from Algorithms.Particle_Swarm.MGALT_PSO_nextGeneration import MGALT_PSO_nextGeneration
from Island_Model.ISL_modSelection import ISL_modSelection
from IPython import get_ipython
from joblib import Parallel, delayed
from Other_Functions.MGALT_bestLegs import MGALT_bestLegs


def MGALT_PSO(bod,const,opts,var,selected,num_mig,num_isl,count_PSO):
    
    ## Pre-Allocation
    bees = np.zeros((opts['PSO']['Npop'],2*opts['PSO']['tspan']))
    dates = np.zeros((opts['PSO']['Npop']*opts['PSO']['tspan'],2*var['transfers']))
    memberDates = np.zeros((1,2*var['transfers']))


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
    f = np.zeros((opts['PSO']['Npop'],1)) 
    f_trans = np.zeros((opts['PSO']['Npop'],var['transfers']))
    avgcost = np.zeros((opts['PSO']['tspan'],1))
    mincost = np.zeros((opts['PSO']['tspan'],1))
    maxcost = np.zeros((opts['PSO']['tspan'],1))
    all_fs = []
    all_fs_trans = np.ndarray((0,var['transfers']))

    
    # Solution Stagnation
    stagnation = 0 
    
    # Number of evals
    nfeval = 0
    
    # Make arrays of 0's for the best solution and cost
    optimal_soln = np.zeros((opts['PSO']['Npop'],np.shape(var['bin'])[0]))
    f_best = 999999999999999999999*np.ones((opts['PSO']['Npop'],1))
    
    ## Generate Initial Population
    
    [popn,bee] = genPopn('PSO',opts,opts['PSO'],var)

    ## Add Shared Solutions if Migration Has Occurred
    
    if num_mig > 0:
       popn = ISL_modReplacement(opts['island'],selected,num_mig-1,num_isl,popn)
       for p in range(opts['PSO']['Npop']):
          bee[p]['pos'][0] = popn[p,:]
    
    ## Iterate Generations
    
    # Progress Display
    get_ipython().magic('clear')
    print('~~~~~~~~~~~ Migrations Completed ',num_mig,'/',opts['island']['Nmig'],'~~~~~~~~~~~')
    print('    Island: ',num_isl+1,'/',opts['island']['Nisl'])
    print('    Algorithm: ',opts['island']['isl_list'][num_isl]) 
    if opts['parallel'] == 'y':
        print('\n    Parallel Processing   \n')
        print('    Watch CPU/RAM usage   \n\n')
    
    for t in range(opts['PSO']['tspan']):
        
        # Calculate Fitness Value for Each Member
        if opts['parallel'] == 'y':
            
            # Display Info

            ll_time = print('\r',' Time Step: ',t+1,'/',opts['PSO']['tspan'],end='') 
            
            # Cost
            fTemp,plot_vars,f_trans = zip(*Parallel(n_jobs=8)(delayed(feval)(popn[member],bod,const,opts,var) for member in range(opts['PSO']['Npop'])))
            fTemp = np.asarray(fTemp)
            
            for member in range(opts['PSO']['Npop']):
                bee[member]['f'] = fTemp[member]
                f[member] = bee[member]['f']
            all_fs = np.append(all_fs,f)
            
            if opts['best_Leg_Tracking'] == 'y':
                all_fs_trans = np.append(all_fs_trans,np.squeeze(f_trans),axis=0)
                var = MGALT_bestLegs(var,popn,np.squeeze(f_trans))

                        
            # Collect results as they become available
            
        else:
        
            # Display Info
            ll_time = print('\r','   Time Step: ',str(t+1),'/',str(opts['PSO']['tspan']) , end = '')
            
            for member in range(opts['PSO']['Npop']):
               bee[member]['f'],pv,f_trans[member,:] = feval(bee[member]['pos'][0],bod,const,opts,var)
               f[member] = bee[member]['f']
            all_fs = np.append(all_fs,f)
            if opts['best_Leg_Tracking'] == 'y':
                all_fs_trans = np.append(all_fs_trans,f_trans,axis=0)
                var = MGALT_bestLegs(var,popn,f_trans)
        
        # See If Fitness Is Best Yet
        for member in range(opts['PSO']['Npop']):
            if f[member] < bee[member]['f_p']:
                bee[member]['p'][0] = bee[member]['pos'][0]
                bee[member]['f_p'] = int(f[member])
        
        # Generate Statistics for Current Time
        avgcost[t] = np.mean(f)
        mincost[t] = min(f)
        maxcost[t] = max(f)
        
        # Organize Best So Far Solutions
        [f_best,optimal_soln,stagnation] = ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation)
            
        # Get dates to see clusters            
        for it in range(np.shape(popn)[0]):
            memberDates[0,:] = popn[it,I]       
            dates[(t)*opts['PSO']['Npop'] + it,:] = 1*memberDates
            
        # Population Members Move, Communicate, & Adjust Velocities
        bee = MGALT_PSO_nextGeneration(bod,opts,opts['PSO'],var,bee)
        
        for p in range(opts['PSO']['Npop']):
            popn[p,:] = bee[p]['pos'][0]
        
        # popn_temp = MGALT_constraintPruning(bod,opts,const,var,popn)
        # for p in range(opts['PSO']['Npop']): 
        #     if popn_temp[p,2] != popn[p,2]:
        #         popn[p,:] = popn_temp[p,:]
        #         bee[p]['pos'][0] = popn_temp[p,:]
        nfeval = nfeval + opts['PSO']['Npop']
    
    
    ## Select Solutions For Sharing
    
    [selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best,optimal_soln,selected,num_isl,num_mig,opts['island'])
    
    ## Eval Info
    
    eval_info = {}
    eval_info['optimal_soln'] = sorted_optimal_soln
    eval_info['f_best'] = sorted_f_best
    eval_info['iterations'] = opts['PSO']['tspan']
    eval_info['maxcost'] = maxcost 
    eval_info['mincost'] = mincost 
    eval_info['avgcost'] = avgcost
    eval_info['total_evals'] = nfeval
    
    return eval_info,selected,all_fs,bees,all_fs_trans
