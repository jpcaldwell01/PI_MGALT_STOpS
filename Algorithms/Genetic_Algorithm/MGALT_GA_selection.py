# FORM: [members,fit] = MGALT_GA_selection(popn,f,opt_GA)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function selects the new members for the Genetic Algorithm
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -popn               (Npop,Nvar)	[float]         [unitless]
# |         Current members of the population
# |     -f                  (Npop,1)  	[float]         [unitless]
# |         The respective cost of each member in 'popn'
# |     -OPT_algo           (1,1)       [struct]        [unitless]
# |         GA option parameters. For a full explination of these 
# |         parameters, see 
# |         "Algorithms/Algorithm_Parameters/parametersGA.m"
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -members        	(Nkeep,Nvar)[float]         [unitless]
# |         Members selected to help produce the next generation
# |     -fit                (Nkeep,1)   [float]         [unitless]
# |         The respective cost of each member in 'members'
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------
import numpy as np
from Other_Functions.randomNum import randomNum

def MGALT_GA_selection(popn,f,opts_algo):

    ## Setup
    [Npop,Nvar] = popn.shape
    
    # Sensible Names to Variables
    elite = opts_algo['elite']          # Number of elite solutions to automatically enter next generation
    N_keep = opts_algo['N_keep']        # Number to select for mating
    T = opts_algo['T']                  # Number of members in each tournament
    method = opts_algo['gen_method']    # Selection method
    threshold = opts_algo['threshold']  # Threshold cost
    
    # Preallocate
    members = np.zeros((elite,np.shape(popn)[1]))
    fit = np.zeros((elite,1))

    # Sort Population from Best --> Worst
    sorted_popn = np.zeros((np.shape(popn)))   # Pre-Allocation
    sorted_f = np.sort(f,axis=0)
    I = np.argsort(f,axis=0)
    
    for v in range(Nvar):
      sorted_popn[:,v] = popn[I,v].flatten()

    ## Selection
    
    if method == 'natural_selection':      # Natural Selection
        '''
        # Elite Members
        if elite
            members(1:elite,:) = sorted_popn(1:elite,:)
            fit(1:elite) = sorted_f(1:elite)
        end

        # Best N Members
        members(elite+1:N_keep+elite,:) = sorted_popn(elite+1:N_keep+elite,:)
        fit(elite+1:N_keep+elite) = sorted_f(elite+1:N_keep+elite)
        '''
    elif method == 'tournament':             # Tournament
        
    	# Which solutions have not been picked yet
        available = list(range(len(sorted_f)))
   
        # Keep "Elite" Members
        if elite:
            for i in range(elite):
                members[i] = sorted_popn[i]
                fit[i] = sorted_f[i]
                available = np.delete(available,i)
   
        # Perform Tournament with Remaining members
        sel = elite + 1
        while len(available) > 1:
            # Select Gladiators
            if len(available) < T:
                competitors = len(available)
            else:
                competitors = T
            
            gladiator = np.zeros((competitors,Nvar))
            gladiator_f = np.zeros((1,competitors))
            
            for i in range(competitors):
                ind = randomNum(1,len(available)-1,'int')
                gladiator[i] = sorted_popn[available[ind]]
                gladiator_f[0,i] = sorted_f[available[ind]]
                available = np.delete(available, ind)
       
            # Fight to the Death
            ind_victor = np.argmin(gladiator_f)
            victor = gladiator[ind_victor]
            victor_f = gladiator_f[0,ind_victor]
       
            # Assign Vicotr to Selected Population
            members = np.append(members,victor.T[None,:],axis=0)
            fit = np.append(fit,victor_f[None,None],axis=0)
            sel = sel + 1
        
    elif method == 'thresholding':           # Thresholding
        '''
    	# Keep "Elite" Members
        if elite
            for i = 1:elite
                members(i,:) = sorted_popn(i,:)
                fit(i) = sorted_f(i)
            end
        end
   
        # Find Solutions that Meet the Threshold
        i = elite + 1
        while sorted_f(i) <= threshold
            members(i,:) = sorted_popn(i,:)
            fit(i) = sorted_f(i)
            i = i + 1
            if i > Npop
                break
            end
        end
        '''
    elif method == 'weighted_random':        # Roulette Wheel
        '''
     	# Keep "Elite" Members
        if elite
            for i = 1:elite
                members(i,:) = sorted_popn(i,:)
                fit(i) = sorted_f(i)
            end
        end
   
        # Pre-Allocate Probability Array
        P = zeros(Npop-elite,1)
    
        # Assign Probabilities for Ranked Weighted Random
        if strcmpi(OPT_algo.weight,'rank')
            den = sum(1:Npop-elite)
            for member = 1:Npop-elite
                P(member) = (Npop-elite-member+1)/den
            end

        # Assign Probabilities for Cost Weighted Random
        elseif strcmpi(OPT_algo.weight,'cost')
            C = -sorted_f(elite+1:end) + max(sorted_f)
            den = sum(C)
            if den ~= 0
                for member = 1:Npop-elite
                    P(member) = abs(C(member)/den)
                end
            else
                P = 1/(Npop-elite)
            end
        end
   
        for pick = 1:N_keep
            # "Spin" The Wheel
            spin = randomNum( 0,sum(P),'dec' )

            # March Through Probability Matrix Until The "Spun" Slice Is Hit
            for i = 1:length(P)
                if spin <= sum(P(1:i))
                    index = i
                    members(pick+elite,:) = sorted_popn(i+elite,:) 
                    fit(pick+elite) = sorted_f(i+elite)
                    break
                end
            end
       
        P(index) = []
        end
        '''
    else:
       
        raise Exception('Incorrect generation method selected. Check the documentation for the specific algorithm to see selection choices.')
    
    return members,fit

