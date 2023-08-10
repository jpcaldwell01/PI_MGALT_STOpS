# FORM: [new_popn,new_f,nfeval] = MGALT_DE_nextGeneration(BOD,CONST,OPT,...
#       OPT_algo,VAR,old_popn,old_f,nfeval,count_DE)
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function creates the next generation of DE popn from the 
# |     old member population and slight perturbations to the old members. 
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
# |         DE option parameters. For a full explination of these 
# |         parameters, see 
# |         "Algorithms/Algorithm_Parameters/parametersDE.m"
# |     -VAR                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |     -old_popn         	(Npop,Nvar)	[float]         [unitless]
# |         The old population members which will be used to create the 
# |         new population members
# |     -old_f              (Npop,1)  	[float]         [unitless]
# |         The respective cost of each member in 'old_popn'
# |     -nfeval           	(1,1)       [int]           [unitless]
# |         The current number of iterations for this DE object
# |     -count_DE           (1,1)       [int]           [unitless]
# |         The current DE island number
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -new_popn           (Npop,Nvar)	[float]         [unitless]
# |         The new population members which generated from the old 
# |         population members and random perturbations
# |     -new_f              (Npop,1)    [float]         [unitless]
# |         The respective cost of each member in 'new_popn'
# |     -nfeval           	(1,1)       [int]           [unitless]
# |         The current number of iterations for this DE object
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np
from Other_Functions.randomNum import randomNum
from Solvers.Indirect import MGALT_IN_FBSM_2D
from Algorithms.Differential_Evolution.MGALT_DE_Target import MGALT_DE_Target
from joblib import Parallel, delayed
from Solvers.Direct import MGALT_DIR_FBSM_2D

def MGALT_DE_nextGeneration(bod,const,opts,opts_algo,var,old_popn,old_f,nfeval,count_DE):

    ## Setup
    
    [Npop,Nvar] = old_popn.shape
    
    if opts['solver'] == 'MGALT_IN_FBSM_2D':
        feval = getattr(MGALT_IN_FBSM_2D,'MGALT_IN_FBSM_2D')
    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
        feval = getattr(MGALT_DIR_FBSM_2D, 'MGALT_DIR_FBSM_2D')
    else:
        raise Exception("Invalid Solver Selection")
    
    # Pre-Allocation
    trial_popn = 1*old_popn
    new_popn = old_popn*0
    trial_f = 1*old_f 
    trial_fs_trans = np.zeros((np.shape(old_f)[0],var['transfers']))
    new_f = old_f*0
    
    ## Base Vector Selection
    
    base_choice = np.zeros((1,Npop),dtype=int)

    if opts_algo['sel_method'] == 'random':             # Random
          '''  
    	# Select Base Vectors So That r0 ~= i
        available = 1:1:Npop
        for i = 1:Npop
            choice = randomNum(1,length(available),'int')
            base_choice(i) = available(choice)
            available(choice) = []
        end
        
        # If Any r0 == i, Fix It
        ctr = 0
        
        while any(base_choice == 1:Npop)
            I = find(base_choice == 1:Npop)
            
            for i = 1:length(I)
                temp = base_choice(I(i))
                
                if I(i) > 1
                    base_choice(I(i)) = base_choice(I(i)-1)
                    base_choice(I(i)-1) = temp
                else
                    base_choice(I(i)) = base_choice(I(i)+1)
                   	base_choice(I(i)+1) = temp
                end
                
            end
            
            ctr = ctr + 1
            
            if ctr > 11
                available = 1:1:Npop
                
                for i = 1:Npop
                    choice = randomNum(1,length(available),'int')
                    base_choice(i) = available(choice)
                    available(choice) = []
                end
                
            end
            
        end
        
        base_vec = old_popn(base_choice,:)
        '''
    elif opts_algo['sel_method'] == 'best_so_far':        # Best so Far
        '''
        [~, ind] = min(old_f)
        base_choice = base_choice + ind
        base_vec = old_popn(base_choice,:)
        '''
    elif opts_algo['sel_method'] == 'random_best_blend':  #Random Best Blend
        
        base_vec = old_popn*0
        ind = np.argmin(old_f)
        best_so_far = old_popn[ind,:]
        
        # Select Base Vectors So That r0 ~= i
        available = list(range(Npop))
        
        for i in range(Npop):
            choice = randomNum(0,len(available)-1,'int')
            base_choice[0,i] = available[choice]
            available = np.delete(available,choice)
        
        # If Any r0 == i, Fix It
        ctr = 0
        
        while np.any(base_choice == list(range(Npop))):
            for i in range(Npop):
                if base_choice[0,i] == i:
                    # np.argwhere(base_choice == list(range(Npop)))
            
                    temp = base_choice[0,i]
                    
                    if i > 1:
                        base_choice[0,i] = base_choice[0,i-1]
                        base_choice[0,i-1] = temp
                    else:
                        base_choice[0,i] = base_choice[0,i+1]
                        base_choice[0,i+1] = temp
            
            ctr = ctr + 1
            
            if ctr > 11:
                available = list(range(Npop))
                
                for i in range(Npop):
                    choice = randomNum(1,len(available),'int')
                    base_choice[i] = available[choice]
                    available = np.delete(available,choice)
            
                ctr = 0
            
        random_choice = old_popn[tuple(base_choice)]
   
        for i in range(Npop):
            base_vec[i] = random_choice[i] + float(np.random.random())*(best_so_far - random_choice[i])
        
    else:
        
        raise Exception('Incorrect selection method selected. Check the documentation for the specific algorithm to see selection choices.')

    ## Generate All Trial Vectors
    
    # Run Through All Target Vectors
    if opts['parallel'] == 'y':
        
        for tv in range(Npop):
            # Generate Target Vectors
            trial_popn[tv] = MGALT_DE_Target(bod,opts,opts_algo,var,old_popn,base_choice,base_vec,Npop,Nvar,tv)
    
        # Run trial solution
        trial_f,plot_vars,trial_fs_trans = zip(*Parallel(n_jobs=8)(delayed(feval)(trial_popn[tv],bod,const,opts,var) for tv in range(opts['DE']['Npop'])))
        trial_f = np.reshape(np.asarray(trial_f),(len(trial_f),1))
        # Collect results as they become available
            
        
    else:
    
        for tv in range(Npop):
    
            # Generate Target Vectors
            trial_popn[tv] = MGALT_DE_Target(bod,opts,opts_algo,var,old_popn,base_choice,base_vec,Npop,Nvar,tv)
    
            # Run trial solution
            trial_f[tv],plot_vars,trial_fs_trans[tv] = feval(trial_popn[tv],bod,const,opts,var)
    
    nfeval = nfeval + Npop
    all_popn = np.vstack((old_popn,trial_popn))
    all_f = np.vstack((old_f,trial_f))
     
    ## Determine Which Members Survive
    
    sorted_popn = all_popn*0 # Pre-Allocation
    
    if opts_algo['surv_method'] == 'natural_selection':  # Natural Selection
        '''
        [sorted_f,I] = sort(all_f)
        
        for v = 1:Nvar
            sorted_popn(:,v) = all_popn(I,v)
        end

        new_popn = sorted_popn(1:Npop,:)
        new_f= sorted_f(1:Npop)
        '''
    elif opts_algo['surv_method'] == 'tournament':         # Tournament
        
        wins = all_f*0 # Pre-Allocation
        
        for i in range(2*Npop):
            available = list(range(0,i)) + list(range(i+1,2*Npop))
            tally = 0
            
            for comp in range(opts_algo['T']):
                c = randomNum(0,len(available)-1,'int')
                competitor = available[c]
                
                if all_f[competitor] <= all_f[i]:
                    tally = tally + 1
                
                available = np.delete(available,c)
            
            wins[i] = tally
        
        I = np.argsort(wins[:,0])
        
        for v in range(Nvar):
            sorted_popn[:,v] = all_popn[I,v]
        
        sorted_f = all_f[I]

        new_popn = sorted_popn[:Npop,:]
        new_f    = sorted_f[:Npop] 
        
    elif opts_algo['surv_method'] == 'weighted_random':    # Weighted Random
        '''
        # Pre-Allocate Probability Array
        P = zeros(2*Npop,1)
        
        # Probability Assignment
        switch OPT_algo.weight
            
            case {'rank'}
                
                den = sum(1:(2*Npop))
            
                for member = 1:(2*Npop)
                    P(member) = ((2*Npop)-member+1)/den
                end
                
            case {'cost'}
                
                normalizer = min(all_f)
            
                for member = 1:(2*Npop)
                    C(member) = all_f(member) - normalizer ##ok<AGROW>
                end

                den = sum(C)

                if den ~= 0
                    for member = 1:(2*Npop)
                        P(member) = abs(C(member)/den)
                    end
                else
                    P = P + 1/(2*Npop)
                end
                
            otherwise
            
                errorPathDisplay()
                fprintf(2,'Incorrect difference vector selected.\n')
                disp('Incorrect difference vector selected.')
                return
            
        end
        
        available = (1:1:(2*Npop))'
        
        for i = 1:Npop
            
            # "Spin" The Wheel
            spin = randomNum( 0,sum(P),'dec' )
            
            # March Through Probability Matrix Until The "Spun" Slice Is Hit
            for j = 1:length(P)
                
                if spin <= sum(P(1:j))
                    new_popn(i,:) = all_popn(available(j),:) 
                    new_f   (i)   = all_f   (available(j))   
                    P(j) = [] 
                    available(j) = []
                    break
                end
                
            end
            
        end
        '''  
    else:
        
        raise Exception('Incorrect survival method selected.')
    
    return new_popn,new_f,nfeval,trial_fs_trans

