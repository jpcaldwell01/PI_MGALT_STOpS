# FORM: [bee] = MGALT_PSO_nextGeneration(BOD,OPT,OPT_algo,VAR,bee)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function creates the next generation of PSO popn from the 
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
# |     -OPT                (1,1)       [struct]        [unitless]
# |         A struct containing constants user options. Contains the save 
# |         folder, ToF values, and more structs containing informaiton 
# |         for the island model, cost parameters, weighting parameters, 
# |         and all of the islands used in the optimization process
# |     -OPT_algo           (1,1)       [struct]        [unitless]
# |         PSO option parameters. For a full explination of these 
# |         parameters, see 
# |         "Algorithms/Algorithm_Parameters/parametersPSO.m"
# |     -VAR                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |     -bee                (Npop,1)	[struct]       	[unitless]
# |         Contains cost, speed position, ect of every bee in the 
# |         population for the old generation
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -bee                (Npop,1)	[struct]       	[unitless]
# |         Contains cost, speed position, ect of every bee in the 
# |         population for the new generation
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np
from Other_Functions.randomNum import randomNum
from Other_Functions.MGALT_fixToF import MGALT_fixToF

def MGALT_PSO_nextGeneration(bod,opts,opt_algo,var,bee):

    ## Setup
    
    # Assign Sensible Names to Variables
    high = np.zeros(np.shape(var['high']))
    low = np.zeros(np.shape(var['low']))

    for it in range(len(var['bin'])):
        high[it] = np.min([var['high'][it],var['highC'][it]])
        low[it] = np.max([var['low'][it],var['lowC'][it]])
        
    Vmax = opt_algo['vmax']*(high-low)
    K = opt_algo['K']
    c1 = opt_algo['cl']
    cmax = opt_algo['cmax']
    Npop = len(bee)
    Nvar = len(bee[0]['pos'][0])
    f_inform = np.zeros((K,1))
    ps = var['ps']
    Index = np.argsort(bee[0]['pos'][0])[::-1] # Shortcut to get indexes of JD regardless of method


    ## Propagate Bees
    
    for b in range(Npop):
       informant = np.zeros((K,1))
       
       # Determine Informants
       for I in range(K):
           while informant[I,0] == 0 or informant[I,0] == b or informant[I,0] == any(informant[:(I-1),0]):
               informant[I,0] = randomNum(0,Npop-1,'int')
               
           f_inform[I,0] = bee[informant[I,0]]['f']
       
       # Determine Best Informant
       ind = np.argmin(f_inform)
       if f_inform[ind,0] < bee[b]['f_g']:
           bee[b]['g'][0] = bee[informant[ind,0]]['pos'][0]
           bee[b]['f_g'] = f_inform[ind,0]
       
       # Move Member & Re-Calculate Velocity
       bee[b]['vel'][0] = c1*bee[b]['vel'][0] + np.random.rand()*cmax*(bee[b]['p'][0]-bee[b]['pos'][0]) + np.random.rand()*cmax*(bee[b]['g'][0]-bee[b]['pos'][0])
       bee[b]['pos'][0] = bee[b]['pos'][0] + bee[b]['vel'][0]
       
       ## Substitution of best legs
       if opts['best_Leg_Tracking'] == 'y':

                # Transfer 1 -> n-1
           for i in reversed(range(1,var['transfers'])):
               if np.random.rand(1) <= ps:
                   bee[b]['pos'][0][Index[i]:Index[i-1]] = var['bestLegMembers']['leg%i' %(var['transfers']-1-i)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1-i)])[0]-1,'int'),Index[i]:Index[i-1]]
                # Transfer n
           if np.random.rand(1) <= ps:
                bee[b]['pos'][0][Index[0]:] = var['bestLegMembers']['leg%i' % (var['transfers']-1)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1)])[0]-1,'int'),Index[0]:]

       # Ensure Bounds Not Exceeded
       for i in range(Nvar):
           
          if bee[b]['vel'][0,i] >  Vmax[i]:
              bee[b]['vel'][0,i] =  Vmax[i]
          
          if bee[b]['vel'][0,i] < -Vmax[i]:
              bee[b]['vel'][0,i] = -Vmax[i]
          
          if bee[b]['pos'][0,i] >  high[i]:
              bee[b]['pos'][0,i] = high[i] 
              bee[b]['vel'][0,i] = -bee[b]['vel'][0,i]
          
          if bee[b]['pos'][0,i] <  low[i]:
              bee[b]['pos'][0,i] = low[i] 
              bee[b]['vel'][0,i] = -bee[b]['vel'][0,i]
          
          if var['bin'][i]:
              bee[b]['pos'][0,i] = round(bee[b]['pos'][0,i])
       
       # Ensure ToF is not broken
       bee[b]['pos'][0] = MGALT_fixToF(bod,opts,var,bee[b]['pos'][0])
    return bee