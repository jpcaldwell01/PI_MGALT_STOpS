# FORM: [trial_popn] = MGALT_DE_Target(BOD,OPT,OPT_algo,VAR,old_popn,...
#       base_choice,base_vec,Npop,Nvar,tv)
# |-----------------------------------------------------------------------
# | NOTES:
# |     -This function creates the target vector for DE_nextGeneration
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
# |     -base_choice     	(1,1)       [float]         [unitless]
# |         Index number for the base choice
# |     -base_vec           (1,1)       [float]         [unitless]
# |         Index number for the base vector selection
# |     -Npop            	(1,1)       [int]           [unitless]
# |         The size of the population
# |     -Nvar             	(1,1)       [int]           [unitless]
# |         The number of variables per member
# |     -tv                 (1,1)       [int]           [unitless]
# |         The current index in "Npop"
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -trial_popn     (1,old_popn)    [float]         [unitless]
# |         The trial population
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

from Other_Functions.randomNum import randomNum
import numpy as np
from Other_Functions.MGALT_fixToF import MGALT_fixToF

def MGALT_DE_Target(bod,opts,opts_algo,var,old_popn,base_choice,base_vec,Npop,Nvar,tv):

    ## Setup
    
    pc = opts_algo['pc']   	# Probability of Crossover
    ps = var['ps']          # Probability of Substitution
    F = opts_algo['F']      	# Scaling Factor
    difference_vec = np.zeros((1,Nvar))
    # Choose Target Vector
    target_vec = old_popn[tv,:]
    high = np.zeros(np.shape(var['high']))
    low = np.zeros(np.shape(var['low']))

    for it in range(len(var['bin'])):
        high[it] = np.min([var['high'][it],var['highC'][it]])
        low[it] = np.max([var['low'][it],var['lowC'][it]])
    
    
    ## Members

    # Choose Two Random Population Members So That r1 ~= r2 ~= i
    choice = randomNum(0,Npop-1,'int')
    while (choice == tv) or (choice == base_choice[0,tv]):
        choice = randomNum(0,Npop-1,'int')
    
    choice1 = choice
    
    choice = randomNum(0,Npop-1,'int')
    while (choice == tv) or (choice == base_choice[0,tv]) or (choice == choice1):
        choice = randomNum(0,Npop-1,'int')
    choice2 = choice
    
    # Compute Weighted Difference Vector
    if opts_algo['F_method'] == 'constant':
        '''
        difference_vec = F*(old_popn[choice1,:] - old_popn[choice2,:])
        '''
    elif opts_algo['F_method'] == 'jitter':

        for i in range(Nvar):
            difference_vec[0,i] = randomNum(F[0],F[1],'dec')*(old_popn[choice1,i] - old_popn[choice2,i])

    elif opts_algo['F_method'] == 'dither':
        '''
        difference_vec = randomNum(F(1),F(2),'dec')*(old_popn(choice1,:) - old_popn(choice2,:))
        '''
    else:

        raise Exception('Incorrect difference vector selected.')
    
    # Add to Base Vector
    mutant_vec = base_vec[tv,:] + difference_vec
    
    trial_popn = 1*target_vec
            
    for i in range(Nvar):
        # Crossover
        if np.random.rand(1) <= pc:
            trial_popn[i] = mutant_vec[0,i]
    
    # Substitution of best legs
    if opts['best_Leg_Tracking'] == 'y':

        I = np.argsort(np.squeeze(mutant_vec))[::-1] # Shortcut to get indexes of JD regardless of method
        for i in reversed(range(1,var['transfers'])):
            if np.random.rand(1) <= ps:
                trial_popn[I[i]:I[i-1]] = var['bestLegMembers']['leg%i' %(var['transfers']-1-i)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1-i)])[0]-1,'int'),I[i]:I[i-1]]
        if np.random.rand(1) <= ps:
            trial_popn[I[0]:] = var['bestLegMembers']['leg%i' % (var['transfers']-1)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1)])[0]-1,'int'),I[0]:]
    
    for i in range(Nvar):
        # Ensure bounds not breached
        if trial_popn[i] > high[i]:
            trial_popn[i] = high[i]
        elif trial_popn[i] < low[i]:
            trial_popn[i] = low[i]
        elif var['bin'][i]:
            trial_popn[i] = round(trial_popn[0,i])
    
    # Ensure ToF is not borken
    trial_popn = MGALT_fixToF(bod,opts,var,trial_popn)
    
    return trial_popn
