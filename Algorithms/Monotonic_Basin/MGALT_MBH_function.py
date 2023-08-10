# FORM: [f_best_return,x_star_best_return,is_viable,neval] = ...
#       MGALT_MBH_function(BOD,CONST,OPT,OPT_algo,VAR,x0,num_mig,count_MBH,loop)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function is the actual implementation of the MBH algorithm. 
# |     The algorithm had to be broken into a function because it would 
# |     not be feasible to fit into "MGALT_MBH" due to organization
# |
# |     -A user who would like solutions where ALL of the transfer 
# |     trajectories are feasible needs to chenge the comment under 
# |     section 3.0 and 4.3 from "any(feasible)" to "all(feasible)". The 
# |     rational for having any was explained in the accompanying thesis, 
# |     as the randomness could potentially alows MBH to find a solution 
# |     where both transfers are considered feasible. As on now, it will 
# |     accept a solution where any transfer is feasible
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
# |     -OPT_algo           (1,1)       [struct]        [unitless]
# |         MBH option parameters. For a full explination of these 
# |         parameters, see 
# |         "Algorithms/Algorithm_Parameters/parametersMBH.m"
# |     -VAR                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |     -x0                 (1,Nvar)    [float]         [unitless]
# |         The initial population which will be evaluated by the solver
# |     -num_mig            (1,1)       [int]           [unitless]
# |         The current migration number
# |     -count_MBH          (1,1)       [int]           [unitless]
# |         The current MBH island number
# |     -loop               (1,1)       [int]           [unitless]
# |         The current loop number
# |         This is if the function is solving an initial set of members 
# |         or a secondary set of members generated from MGALT_MBH_cluster
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -f_best             (1,Nvar)    [float]         [unitless]
# |         The respective cost of each member in 'x0'
# |     -member_final       (1,Nvar)    [float]         [unitless]
# |         The new perturbed member consisting of the original member and
# |         the perturbations added to it
# |     -is_viable          (1,1)       [boolean]   	[unitless]
# |         The new perturbed member consisting of the original member and
# |         the perturbations added to it
# |     -nfeval           	(1,1)       [int]           [unitless]
# |         The current number of iterations for this MBH object
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------
from Solvers.Indirect import MGALT_IN_FBSM_2D
from Solvers.Direct import MGALT_DIR_FBSM_2D
import numpy as np
from Algorithms.Monotonic_Basin.MGALT_MBH_isFeasible import MGALT_MBH_isFeasible
from Algorithms.Monotonic_Basin.MGALT_MBH_basin import MGALT_MBH_basin
import math
from joblib import Parallel, delayed
import gc

def MGALT_MBH_function(bod,const,opts,opt_algo,var,x0,num_mig,count_MBH,loop,Nkeep):

    ## Initialize
    
    # Number of iterations
    num_iter_global = len(x0)     # Number of iterations for the outer loop
    
    # Preallocate Cost Parameters
    f = []
    plot_vars = ()    
    all_fs_trans = np.ndarray((0,var['transfers']))
    f_best_return = np.zeros((Nkeep,1))
    x_star_best_return = np.zeros((num_iter_global,np.shape(x0)[1]))
    is_viable = [False for x in range(num_iter_global)]
    
    # Perallocate the plot_vars struct. Need to run 1 execution to get field
    # names
    if opts['solver'] == 'MGALT_IN_FBSM_2D':
        feval = getattr(MGALT_IN_FBSM_2D, 'MGALT_IN_FBSM_2D')
        Ic = [0] + [11*x-1 for x in range(1,var['transfers'])] + [-1] + [11*x for x in range(1,var['transfers'])] # Cluster pruning indexes

    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
        feval = getattr(MGALT_DIR_FBSM_2D, 'MGALT_DIR_FBSM_2D')
        Ic = [0] + [(opts['thrust']['Nseg']*2+5)*x-1 for x in range(1,var['transfers'])] + [-1] + [(opts['thrust']['Nseg']*2+5)*x for x in range(1,var['transfers'])] # Cluster pruning indexes

    else:
        raise Exception("Invalid Solver Selection")
        
    # Vars = feval(x0[0,:],bod,const,opts,var)[1]
    # name = list(Vars.keys())
    # plot_vars = [{}]*num_iter_global
    # for i0 in range(num_iter_global):
    #     plot_vars[i0] = dict.fromkeys(name,[])
    
    ## feval
        
    if opts['parallel'] == 'y':
        
        # Display info
        ll_par = print('\n    Parallel Processing   \n')
        ll_RAM = print('    Watch CPU/RAM usage   \n\n')
    
        # Chunks to break parallel into to prevent a lot of RAM usage at once
        chunk_nums = 500      # Set by Malloy, 1000 uses ~2.0GB of RAM per pass in MATLAB
        chunk_whole = math.floor(num_iter_global/chunk_nums)
        chunk_rems = num_iter_global % chunk_nums
        
        # Run through all searches minus mod
        for chunk_iter in range(chunk_whole):
            
            # Preallocate temp vars
            # f_temp = np.zeros((chunk_nums,1))
            # for i0 in range(len(name)):
            #     plot_vars_temp[chunk_nums][name[i0]] = []
            
            # Display info
            ll_loop1 = print('\r','    Search: ',chunk_iter*chunk_nums+1,'-',(chunk_iter+1)*chunk_nums,' /',num_iter_global,end='')
            # ll_loop2 = print('\r','    Search: %1.0f\n',end='' % chunk_iter*chunk_nums+2)
            
            # Perform parfeval on chunk n
            # data_temp[chunk_num] = feval(OPT.ppool,fh,2,x0(((chunk_iter-1)*chunk_nums)+chunk_num,:),BOD,CONST,OPT,VAR)
            f_temp,plot_vars_temp,f_trans = zip(*Parallel(n_jobs=8)(delayed(feval)(x0[(bool(chunk_iter)*Nkeep)+chunk_num,:],bod,const,opts,var) for chunk_num in range(chunk_nums)))
            f = np.append(f,np.asarray(f_temp))
            plot_vars = np.append(plot_vars,plot_vars_temp)
            if opts['best_Leg_Tracking'] == 'y':
                all_fs_trans = np.append(all_fs_trans,np.squeeze(f_trans),axis=0)

            
            # Clear Memory
            I = np.argsort(f)
            f = np.delete(f,I[Nkeep:])
            plot_vars = np.delete(plot_vars,I[Nkeep:])
            if opts['best_Leg_Tracking'] == 'y':
                all_fs_trans = np.delete(all_fs_trans,I[Nkeep:],0)
            x0 = np.delete(x0,I[Nkeep:],0)
            
            del f_temp,plot_vars_temp,f_trans
            gc.collect()
        
        # Run through all remaining values
        if bool(chunk_rems):
                
            # Display info
            ll_loop1 = print('\r','    Search: ',chunk_whole*chunk_nums+1,'-',chunk_whole*chunk_nums+chunk_rems,' / ',num_iter_global,end='')
    
            # Perform parfeval on remainder
            f_temp,plot_vars_temp,f_trans = zip(*Parallel(n_jobs=8)(delayed(feval)(x0[Nkeep+chunk_rem,:],bod,const,opts,var) for chunk_rem in range(chunk_rems)))
            f = np.append(f,np.asarray(f_temp))
            plot_vars = np.append(plot_vars,plot_vars_temp)
            if opts['best_Leg_Tracking'] == 'y':
                all_fs_trans = np.append(all_fs_trans,np.squeeze(f_trans),axis=0)

            # Clear Memory
            I = np.argsort(f)
            f = np.delete(f,I[Nkeep:])
            plot_vars = np.delete(plot_vars,I[Nkeep:])
            if opts['best_Leg_Tracking'] == 'y':
                all_fs_trans = np.delete(all_fs_trans,I[Nkeep:],0)
            x0 = np.delete(x0,I[Nkeep:],0)
            del f_temp,plot_vars_temp,f_trans
        
    else:
        
        # Display info        
        # Cost
        i1 = 0
        count = 1
        while i1 < len(x0[:,0])-1:
    
            # Search Display
            ll_loop = print('\r','    Search: ',count,'/',num_iter_global,end='')
            count = count+1
            
            # ***2.0 1st run of NLP***
            #Run the NLP problem solver to find x_star using the initial guess x0
            [f_temp,plot_vars_temp,f_trans] = feval(x0[i1,:],bod,const,opts,var)
            f = np.append(f,np.asarray(f_temp))
            plot_vars = np.append(plot_vars,plot_vars_temp)
            all_fs_trans = np.append(all_fs_trans,f_trans,axis=0)
            
            i1 = i1+1
            
            # Clear Memory
            if len(f) > Nkeep:
                I = np.argsort(f)
                f = np.delete(f,I[Nkeep:])
                plot_vars = np.delete(plot_vars,I[Nkeep:])
                all_fs_trans = np.delete(all_fs_trans,I[Nkeep:],0)
                x0 = np.delete(x0,I[Nkeep:],0)
                i1 = i1-1
    
    neval = num_iter_global
    
    ## Check Feasibility
    
    ll_feas = print('\n ---Checking Feasibility---\n')
    
    per_feas = opt_algo['per_feas'][num_mig]	# Feasibility percentage for SOI
    
    search_index = []
    count_index = 0
    
    for i2 in range(len(f)):
        
        # ***3.0 Check x_star feasibility***
        feasible = MGALT_MBH_isFeasible(bod,const,opts,opt_algo,var,plot_vars[i2],per_feas)[0]
        
        if opt_algo['feas_check'](feasible):
            
            search_index[count_index] = i2
            count_index = count_index+1
            
            is_viable[i2] = True
    
        else:
            
            # No basin was found, try the next member
            f_best_return[i2] = f[i2]
            x_star_best_return[i2,:] = x0[i2,:]
            is_viable[i2] = False
    
    # If nothing feasible was found, return
    if len(search_index) == 0:
        return f_best_return,x_star_best_return,is_viable,neval,Ic,all_fs_trans

    ## If Feasibility, Check Basins
    
    if loop == 1:
        disp_info = '\n ---Potential Basin Found---'
    elif loop == 2:
        disp_info = '\n ---Exploring Minima---'
    
    # A basin has been found, start the iteration process to explore it
    ll_MBH = print(disp_info)
    
    if opts['parallel'] == 'y':
        '''
        # Preallocate
        fh = @MGALT_MBH_basin
        b_data(1:size(search_index,2)) = parallel.FevalFuture
        
        # Display info
        ll_par = fprintf('\n    Parallel Processing   \n')
        ll_RAM = fprintf(2,'    Watch CPU/RAM usage   \n\n')
        
        # feval
        for i3 = 1:size(search_index,2)
            b_data(i3) = parfeval(OPT.ppool,fh,2,BOD,CONST,OPT,OPT_algo,VAR,...
                f(search_index(i3)),x0(search_index(i3),:),num_mig,count_MBH,loop)
        end
        
         # Collect results as they become available
        for i4 = 1:size(search_index,2)
            [p_index,p_f,p_x] = fetchNext(b_data)
            p_f_best(p_index) = p_f
            p_x_best(p_index,:) = p_x
        end
        
        # Put the data into the correct index for func return
        for i5 = 1:size(search_index,2)
            f_best_return(search_index(i5)) = p_f_best(i5)
            x_star_best_return(search_index(i5),:) = p_x_best(i5,:)
        end
    
        # Clear display
        fprintf(repmat('\b',1,ll_RAM))
        fprintf(repmat('\b',1,ll_par))
       ''' 
    else:
        
        for i3 in range(len(search_index)):
    
            [f_best_return[search_index[i3]],x_star_best_return[search_index[i3],:]] = MGALT_MBH_basin(bod,const,opts,opt_algo,var,f[search_index[i3]],x0[search_index[i3],:],num_mig,count_MBH,loop)
                
    # Remove "Potential Basin Found"
    # fprintf(repmat('\b',1,ll_MBH))
    
    return f_best_return,x_star_best_return,is_viable,neval,Ic,all_fs_trans
