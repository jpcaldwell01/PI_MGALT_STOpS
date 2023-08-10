# FORM: [child1,child2] = MGALT_GA_mating(BOD,OPT,opts_algo,VAR,parent1,parent2)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function allows two different parents to mate and produce 
# |     two different offspring. Used in the Genetic Algorithm (GA)
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
# |     -opts_algo           (1,1)       [struct]        [unitless]
# |         GA option parameters. For a full explination of these 
# |         parameters, see 
# |         "Algorithms/Algorithm_Parameters/parametersGA.m"
# |     -VAR                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |     -parent1          	(1,Nvar)	[float]         [unitless]
# |         First parent member of the population
# |     -parent2            (1,Nvar)  	[float]         [unitless]
# |         Second parent member of the population
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -child1             (1,Nvar)	[float]         [unitless]
# |         First child member of the new population
# |     -child2             (1,Nvar)	[float]         [unitless]
# |         Second child member of the new population
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np
from Other_Functions.MGALT_fixToF import MGALT_fixToF
from Other_Functions.randomNum import randomNum

def MGALT_GA_mating(bod,opts,opts_algo,var,parent1,parent2):

    ## Setup
    
    Nvar = len(parent1)  
    child1 = parent1
    child2 = parent2 # In case neither crossover nor mutation occurs
    l = np.zeros((Nvar,1))
    
    ## Mating
    
    # Sensible Names to Variables
    high = np.zeros(np.shape(var['high']))
    low = np.zeros(np.shape(var['low']))

    for it in range(len(var['bin'])):
        high[it] = np.min([var['high'][it],var['highC'][it]])
        low[it] = np.max([var['low'][it],var['lowC'][it]])
    bins = var['bin']
    pc = opts_algo['pc']
    pm = opts_algo['pm']
    ps = var['ps']
    
    if opts_algo['mate_method'] in ['random_crossover','uniform_crossover']:
            
            # For 3 digit availability, slide variables over 3 digits
            fac = 1000
            p_bin1 = ''
            p_bin2 = ''
       
            # If there is a negative number, slide all variables up by that number
            min_num = np.min([parent1,parent2])
            if min_num < 0:
                parent1 = parent1 + np.abs(min_num)
               	parent2 = parent2 + np.abs(min_num)
            
    
            # Get all binary arrays Make matching arrays equal in length
            for i in range(Nvar):
                p1 = np.binary_repr(int(parent1[i]*fac))
                
                if np.isinf(parent1[i]*fac) or np.isinf(parent2[i]*fac):
                    print('u')
                
                
                p2 = np.binary_repr(int(parent2[i]*fac))
                
                while len(p1) < len(p2):
                    p1 = '0'+ p1
                
                while len(p2) < len(p1):
                    p2 = '0'+ p2
                
                l[i] = len(p1) ##ok<AGROW>
                p_bin1 = p_bin1 + p1
                p_bin2 = p_bin2 + p2 ##ok<AGROW>
            
            p_bin1List = 1*list(p_bin1)
            p_bin2List = 1*list(p_bin2)
            kid1 = 1*list(p_bin1)
            kid2 = 1*list(p_bin2)
            
            
            # Perform the bit crossover
            N = len(p_bin1)
            if opts_algo['cross_points'] > (N-1):
                opts_algo['cross_points'] = N-1 
                
            if opts_algo['mate_method'] == 'uniform_crossover':
                opts_algo['cross_points'] = N-1 
                
            available = list(range((N-1)))
            cut = np.zeros((opts_algo['cross_points']+1,1))
            
            for i in range(opts_algo['cross_points']): # Determine crossover cutoff points
                choice = randomNum(1,len(available)-1,'int')
                cut[i] = available[choice] ##ok<AGROW>
                available = np.delete(available,choice)
            
            cut = np.sort(cut,axis=0)
            for i in range(len(cut)-1):
                
                array = [cut[i],cut[i+1]]
                
                if np.random.rand(1) <= pc:
                    kid1[int(array[0]):int(array[-1])] = p_bin2List[int(array[0]):int(array[-1])] 
                    
                if np.random.rand(1) <= pc:
                    kid2[int(array[0]):int(array[-1])] = p_bin1List[int(array[0]):int(array[-1])] 
            
            # Turn the binary strings back into variable arrays
            for i in range(Nvar):
                
                array = [0,int(l[i])-1]
                
                if i != 0:
                    array = [val+np.sum(l[:i]) for val in array]
                
                child1[i] = int(''.join(kid1[int(array[0]):int(array[-1])+1]),2)/fac
                child2[i] = int(''.join(kid2[int(array[0]):int(array[-1])+1]),2)/fac
       
            # If the number were slid to avoid negatives, slide them back
            if min_num < 0:
                child1 = child1 - np.abs(min_num)
                child2 = child2 - np.abs(min_num)
            
    elif opts_algo['mate_method'] == 'blending':
           ''' 
        	child1 = parent1 
            child2 = parent2
            if parent1 ~= parent2 
                for i = 1:Nvar
                    if rand <= pc
                        beta1 = randomNum(-opts_algo.OB,1+opts_algo.OB,'dec')
                       	child1(i) = beta1*parent1(i) + (1-beta1)*parent2(i)
                    end
                    if rand <= pc
                        beta2 = randomNum(-opts_algo.OB,1+opts_algo.OB,'dec')
                       	child2(i) = beta2*parent1(i) + (1-beta2)*parent2(i)
                    end
                end
       
            # If both parents are the same, there is a chance to 'bump' the
            # solution a little bit to see if it improves, instead of just
            # having both children be the same as the parents
            else
                range = high - low
                for i = 1:Nvar
                    if rand <= pc
                        beta1 = randomNum(-opts_algo.OB,opts_algo.OB,'dec')
                      	child1(i) = (beta1*range(i)) + parent1(i)
                    end
                    if rand <= pc
                        beta2 = randomNum(1-opts_algo.OB,1+opts_algo.OB,'dec')
                      	child2(i) = (beta2*range(i)) + parent1(i)
                    end
                end
            end
           ''' 
    else:
            
        raise Exception('Incorrect mating method selected. Check the documentation for the specific algorithm to see selection choices.')
    
    ## Mutation
    
    # Mutate All Variables
    if np.random.rand(1) <= pm:
       for i in range(Nvar):
          child1[i] = (high[i]-low[i])*float(np.random.rand(1)) + low[i]

    if np.random.rand(1) <= pm:
       for i in range(Nvar):
          child2[i] = (high[i]-low[i])*float(np.random.rand(1)) + low[i]
          
    ## Substitution of best legs
    if opts['best_Leg_Tracking'] == 'y':
        I = np.argsort(child1)[::-1] # Shortcut to get indexes of JD regardless of method
        
             # Transfer 1 -> n-1
        for i in reversed(range(1,var['transfers'])):
            if np.random.rand(1) <= ps:
                child1[I[i]:I[i-1]] = var['bestLegMembers']['leg%i' %(var['transfers']-1-i)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1-i)])[0]-1,'int'),I[i]:I[i-1]]
    
        for i in reversed(range(1,var['transfers'])):
            if np.random.rand(1) <= ps:
                child2[I[i]:I[i-1]] = var['bestLegMembers']['leg%i' %(var['transfers']-1-i)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1-i)])[0]-1,'int'),I[i]:I[i-1]]
                
        if np.random.rand(1) <= ps:
            child1[I[0]:] = var['bestLegMembers']['leg%i' % (var['transfers']-1)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1)])[0]-1,'int'),I[0]:]
        if np.random.rand(1) <= ps:
            child2[I[0]:] = var['bestLegMembers']['leg%i' % (var['transfers']-1)][randomNum(0,np.shape(var['bestLegMembers']['leg%i' %(var['transfers']-1)])[0]-1,'int'),I[0]:]

    
    ## Ensure Bounds Not Broken and Binary Stays Binary
    
    for i in range(Nvar):
       if child1[i] > high[i]:
           child1[i] = high[i]
       
       if child1[i] < low[i]:
           child1[i] = low[i]
       
       if child2[i] > high[i]:
           child2[i] = high[i] 
       
       if child2[i] < low[i]:
           child2[i] = low[i]
       
       if bins[i]:
           child1[i] = round(child1[i])
           child2[i] = round(child2[i])
    
    ## Correct for ToF 
    # Added by Malloy to account for the incorrect ToF variables between the
    # members
    
    child1 = MGALT_fixToF(bod,opts,var,child1)
    child2 = MGALT_fixToF(bod,opts,var,child2)
    
    return child1,child2
