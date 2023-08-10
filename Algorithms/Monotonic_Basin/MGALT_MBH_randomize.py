# FORM: [x_prime] = MGALT_MBH_randomize(BOD,OPT,VAR,x_current,per_rand)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function creates random perturbations within x_current and
# |     returns it as x_prime
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -BOD                (1,1)       [struct]        [unitless]
# |         A struct containing information pertaining to the planetary
# |         bodies. Contains list of bodies, launch windows and ToF, and 
# |         planetary R/V/JD vectors. This struct has dynamic fields and 
# |         will adapt to contain only the necesary information
# |     -OPT                (1,1)       [struct]        [unitless]
# |         A struct containing constants user options. Contains the save 
# |         folder, ToF values, and more structs containing informaiton 
# |         for the island model, cost parameters, weighting parameters, 
# |         and all of the islands used in the optimization process
# |     -VAR                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |     -x_current          (1,Nvar)    [float]         [unitless]
# |         The current population which will be perturbed
# |     -per_rand           (1,1)       [float]     	[unitless]
# |         Random percentage number for generating new results
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -x_prime            (1,Nvar)    [float]         [unitless]
# |         The final perturbed population
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np
from Other_Functions.randomNum import randomNum
from Other_Functions.MGALT_fixToF import MGALT_fixToF

def MGALT_MBH_randomize(bod,opts,var,x_current,per_rand):

    ## Perform Randomization
    
    # Perallocate
    x_prime = np.zeros((len(x_current),1))
    
    # Randomization
    if opts['solver'] in ['MGALT_DIR_FSM_2D','MGALT_DIR_FBSM_2D']:
            '''
            # Apply perturbations on every member
            for i1 in range(len(x_prime)):
                random_gen = x_current[i1]*per_rand
                
                if var['bin'][i1]:   # Binary variables
                    x_prime[i1] = randomNum(x_current[i1]-random_gen,x_current[i1]+random_gen, 'int')
                else:
                    x_prime[i1] = randomNum(x_current[i1]-random_gen,x_current[i1]+random_gen, 'dec')
            
            # Ensure the first thrust = 1 or 0
            switch OPT.thrust.thrust_method
                
                case {'constant_thrust','equation_thrust'}
                    
                    x_prime(2) = 1     # Ensure that always thrusting when leaving departure planet
                    
                case {'variable_thrust'}
                    
                    # Acts as a pass
                    
                otherwise
                    
                    errorPathDisplay()
                    fprintf(2,'Incorrect thrust method selected.\n')
                    return
                    
            end
            '''
    elif opts['solver'] in ['MGALT_IN_FSM_2D','MGALT_IN_FBSM_2D']:

        # Apply perturbations on every member
        for i1 in range(len(x_prime)):
            random_gen = x_current[i1]*per_rand
            x_prime[i1] = randomNum(x_current[i1]-random_gen,x_current[i1]+random_gen, 'dec')
        
    else:
        raise Exception("Invalid Solver Selection")
    
    
    # Ensure the bounds of the search space haven't been exceeded
    for i5 in range(len(x_prime)):
        
        if x_prime[i5] < var['low'][i5]:
            x_prime[i5] = var['low'][i5]
        
        if x_prime[i5] > var['high'][i5]:
            x_prime[i5] = var['high'][i5]
    
    # Ensure that ToF is not broken
    x_prime = MGALT_fixToF(bod,opts,var,x_prime)
    
    return x_prime
