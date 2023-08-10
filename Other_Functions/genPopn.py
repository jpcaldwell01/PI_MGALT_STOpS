# FORM: [popn,bee] = genPopn(algorithm,OPT,opt_algo,VAR)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Function designed to randomly generate an initial population for 
# |     one of the algorithms. Each algorithms has a specificway to 
# |     generate the population, which is specified in the switch/case 
# |     by the algorithm string. This funciton also ensures that the Tof 
# |     bounds are not broken from the start, so T1_JD+T1_ToF == T2_JD
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -algorithm          (1,n)       [string]        [unitless]
# |         The solver (cost function) handle which is used to determine 
# |         the correct solver to use. This string is used in an ode45 
# |         function call and is case specific for the function.m name
# |     -OPT                (1,1)       [struct]        [unitless]
# |         A struct containing constants user options. Contains the save 
# |         folder, ToF values, and more structs containing informaiton 
# |         for the island model, cost parameters, weighting parameters, 
# |         and all of the islands used in the optimization process
# |     -OPT_algo           (1,1)       [struct]        [unitless]
# |         Algo option parameters
# |     -VAR                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -popn           (Nshared,Nvar)	[float]         [unitless]
# |         A randomely generated array of members from the initalization 
# |         of an island
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

def genPopn(algorithm,opts,opt_algo,var):

    ## Initials
    
    Nvar = len(var['bin'])       	# Number of Variables
    high = np.zeros(np.shape(var['high']))
    low = np.zeros(np.shape(var['low']))

    for it in range(len(var['bin'])):
        high[it] = np.min([var['high'][it],var['highC'][it]])
        low[it] = np.max([var['low'][it],var['lowC'][it]])
    
    
    ## Pick which solver/algorithm to generate for
    
    if opts['solver'] in ['LT_DIR_FSM_2D','LT_IN_FSM_2D']:
        
        segment = Nvar/var['transfers']	# Number of segments for the population
        
        if algorithm in ['GA','DE']:

            # Initials
            bee = 0

            # Pre-Allocate with Random Normalized Values
            popn = np.random.rand(opt_algo['Npop'],Nvar) #population

            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 in range(opt_algo['Npop']):
                for i2 in range(Nvar):
                    if not var['bin'][i2]:      # Non-binary variables
                        popn[i1,i2] = (var['high'][i2]-var['low'][i2])*popn[i1,i2] + var['low'][i2]  
                    elif var['bin'][i2]:   # Binary variables
                        popn[i1,i2] = np.round(popn[i1,i2])

            # Break the popn into transfer segments
            # ie, popn is 75x44 with 2 transfers, so the first stansfer 
            # segment is popn(:,1:21) and the second segment is popn(:,22:end)
            # This needs to account for the time of flight, as the second
            # population segment can't start before the first one, so the
            # transfer departure JD must be adjusted
            for i3 in range(opt_algo['Npop']):
                for i4 in range(var['transfers']-1):
                    start = popn[i3,(segment*i4)-segment]
                    tof = popn[i3,segment*i4]
                    popn[i3,(segment*i4)+1] = start+tof

        
        elif algorithm == 'PSO':
            '''
            # Initials
            Vmax = opt_algo['vmax']*(var['high']-var['low'])

            # Pre-Allocate
            popn = np.zeros((opt_algo['Npop'],Nvar))
            bee = {'pos': np.zeros((1,Nvar)), 'vel': np.zeros((1,Nvar)), 'f': 0, 'f_p':999999999999999, 'p' : np.zeros((1,Nvar)), 'f_g':999999999999999, 'g': np.zeros((1,Nvar))} 
            print()
            # bee[:opts_algo['Npop']] = {'pos': np.zeros((1,Nvar)), 'vel': np.zeros((1,Nvar)), 'f': 0, 'f_p':999999999999999, 'p' : np.zeros((1,Nvar)), 'f_g':999999999999999, 'g': np.zeros((1,Nvar))} 
            
            # Randomize the Values
            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 = 1:OPT_algo.Npop
                for i2 = 1:Nvar
                    bee(i1).pos(i2) = randomNum(VAR.low(i2),VAR.high(i2),'dec')
                    bee(i1).vel(i2) = randomNum(-Vmax(i2),Vmax(i2),'dec')
                    if VAR.bin(i2)
                        bee(i1).pos(i2) = round(bee(i1).pos(i2))
                    end
                end
            end

            # Break the popn into transfer segments
            # ie, popn is 75x40 with 2 transfers, so the first stansfer 
            # segment is popn(:,1:20) and the second segment is popn(:,21:end)
            for i3 = 1:OPT_algo.Npop
                for i4 = 1:(VAR.transfers-1)
                    start = bee(i3).pos((segment*i4)-segment+1)
                    tof = bee(i3).pos(segment*i4)
                    bee(i3).pos((segment*i4)+1) = start+tof
                end
                popn(i3,:) = bee(i3).pos
            end

        elif algorithm == 'MBH':

            # Initials
            bee = 0

            # Pre-Allocate with Random Normalized Values
            popn = rand(OPT_algo.N1_Outer,Nvar) # Population

            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 = 1:OPT_algo.N1_Outer
                for i2 = 1:Nvar
                    if ~VAR.bin(i2)      # Non-binary variables
                        popn(i1,i2) = (VAR.high(i2)-VAR.low(i2))*popn(i1,i2) + VAR.low(i2) 
                    elseif VAR.bin(i2)   # Binary variables
                        popn(i1,i2) = round(popn(i1,i2))
                    end
                end
            end

            # Break the popn into transfer segments
            # ie, popn is 75x40 with 2 transfers, so the first stansfer 
            # segment is popn(:,1:20) and the second segment is popn(:,21:end)
            for i3 = 1:OPT_algo.N1_Outer
                for i4 = 1:(VAR.transfers-1)
                    start = popn(i3,(segment*i4)-segment+1)
                    tof = popn(i3,segment*i4)
                    popn(i3,(segment*i4)+1) = start+tof
                end
            end

        else:

            raise Exception("Invalid Algorithm Selection")
            return
        '''
    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
        
        if algorithm in ['GA','DE']:
            
            # Initials
            bee = 0

            # Pre-Allocate with Random Normalized Values
            popn = np.random.rand(opt_algo['Npop'],Nvar) #population

            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 in range(opt_algo['Npop']):
                for i2 in range(Nvar):
                    if not var['bin'][i2]:      # Non-binary variables
                        popn[i1,i2] = (high[i2]-low[i2])*popn[i1,i2] + low[i2]  
                    elif var['bin'][i2]:   # Binary variables
                        popn[i1,i2] = np.round(popn[i1,i2])
            
            # Break the popn into transfer segments
            if var['transfers'] > 1:
                
                # How many segments is the transfers broken into
                seg = len(popn.transpose()) - (2 + ((var['transfers']-1)*5))    # Disregarding the misc info, how many thrust/angle
                seg = seg/var['transfers']                             # Thrust/angler per transfer
                seg = int(seg/2)                                            # How many segments
                
                for i3 in range(opt_algo['Npop']):
                    
                    # From Departure to flyby 1
                    start = popn[i3,0]
                    tof = popn[i3,seg*2+4]
                    popn[i3,seg*2+5] = start+tof
                
                    # From flyby 2 to n
                    for i4 in range(2,var['transfers']):
                        start = popn[i3,((i4-1)*((2*seg)+5))]
                        tof = popn[i3,((i4)*((2*seg)+5))-1]
                        popn[i3,(i4)*((2*seg)+5)] = start+tof

         
        elif algorithm in 'PSO':
            
            # Initials
            Vmax = opt_algo['vmax']*(high-low)
            bee = {}
            
            # Pre-Allocate
            popn = np.zeros((opt_algo['Npop'],Nvar))
            # bee(1:OPT_algo.Npop) = struct('pos',zeros(1,Nvar), 'vel',zeros(1,Nvar), 'f',0,'f_p',999999999999999, 'p',zeros(1,Nvar), 'f_g',999999999999999, 'g',zeros(1,Nvar)) 
            for i in range(opt_algo['Npop']):
                bee[i] = {'pos':np.zeros((1,Nvar)), 'vel':np.zeros((1,Nvar)), 'f':0,'f_p':999999999999999, 'p':np.zeros((1,Nvar)), 'f_g':999999999999999, 'g':np.zeros((1,Nvar))}
            
            # Randomize the Values
            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 in range(opt_algo['Npop']):
                for i2 in range(Nvar):
                    bee[i1]['pos'][0,i2] = randomNum(low[i2],high[i2],'dec')
                    bee[i1]['vel'][0,i2] = randomNum(-Vmax[i2],Vmax[i2],'dec')
                    if var['bin'][i2]:
                        bee[i1]['pos'][i2] = round(bee[i1]['pos'][i2])

            for i1 in range(opt_algo['Npop']):
                for i2 in range(Nvar):
                    if var['bin'][i2] == 0:      # Non-binary variables
                        popn[i1,i2] = (high[i2] - low[i2])*popn[i1,i2] + low[i2]
                    elif var['bin'][i2]:   # Binary variables
                        popn[i1,i2] = round(popn[i1,i2])

            # Break the popn into transfer segments
            if var['transfers'] > 1:
                
                # How many segments is the transfers broken into
                seg = np.shape(popn)[1] - (2 + ((var['transfers']-1)*5))    # Disregarding the misc info, how many thrust/angle
                seg = seg/var['transfers']                             # Thrust/angler per transfer
                seg = int(seg/2)                                            # How many segments
                
                for i3 in range(opt_algo['Npop']):
                    
                    # From Departure to flyby 1
                    start = bee[i3]['pos'][0,0]
                    tof = bee[i3]['pos'][0,seg*2+4]
                    bee[i3]['pos'][0,seg*2+5] = start+tof
                
                    # From flyby 2 to n
                    for i4 in range(1,var['transfers']-1):
                        start = bee[i3]['pos'][0,i4*((2*seg)+5)]
                        tof = bee[i3]['pos'][0,(i4+1)*((2*seg)+5)-1]
                        bee[i3]['pos'][0,(i4*((2*seg)+5))] = start+tof
                        
                    popn[i3,:] = bee[i3]['pos']
        
        elif algorithm == 'MBH':
            
            # Initials
            bee = 0

            # Pre-Allocate with Random Normalized Values
            popn = np.random.rand(opt_algo['N1_Outer'],Nvar) #population

            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 in range(opt_algo['N1_Outer']):
                for i2 in range(Nvar):
                    if var['bin'][i2] == 0:      # Non-binary variables
                        popn[i1,i2] = (high[i2] - low[i2])*popn[i1,i2] + low[i2]  
                    elif var['bin'][i2]:   # Binary variables
                        popn[i1,i2] = round(popn[i1,i2])
            
            # Break the popn into transfer segments
            if var['transfers'] > 1:
                
                # How many segments is the transfers broken into
                seg = np.shape(popn)[1] - (2 + ((var['transfers']-1)*5))    # Disregarding the misc info, how many thrust/angle
                seg = seg/var['transfers']                             # Thrust/angler per transfer
                seg = int(seg/2)                                            # How many segments
                
                for i3 in range(opt_algo['N1_Outer']):
                    
                    # From Departure to flyby 1
                    start = popn[i3,0]
                    tof = popn[i3,seg*2+4]
                    popn[i3,seg*2+5] = start+tof
                
                    # From flyby 2 to n
                    for i4 in range(1,var['transfers']-1):
                        start = popn[i3,i4*((2*seg)+5)]
                        tof = popn[i3,(i4+1)*((2*seg)+5)-1]
                        popn[i3,(i4+1)*((2*30)+5)] = start+tof
            
        else:
            
            raise Exception("Invalid Island Selection")
            
        
    elif opts['solver'] == 'MGALT_IN_FBSM_2D':
            
        if algorithm in ['GA','DE']:
            
            # Initials
            bee = 0

            # Pre-Allocate with Random Normalized Values
            popn = np.random.rand(opt_algo['Npop'],Nvar) #population
            
            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 in range(opt_algo['Npop']):
                for i2 in range(Nvar):
                    if not var['bin'][i2]:      # Non-binary variables
                        popn[i1,i2] = (high[i2]-low[i2])*popn[i1,i2] + low[i2]  
                    elif var['bin'][i2]:   # Binary variables
                        popn[i1,i2] = np.round(popn[i1,i2])
                    
            
            # Break the popn into transfer segments
            for i3 in range(opt_algo['Npop']):
                for i4 in range(var['transfers']-1):
                    start = popn[i3,(i4+1)*11-11]
                    tof = popn[i3,(i4+1)*11-1]
                    popn[i3,(i4+1)*11] = start+tof
                
            
            
        elif algorithm == 'PSO':
            
            # Initials
            Vmax = opt_algo['vmax']*(high-low)
            bee = {}
            
            # Pre-Allocate
            popn = np.zeros((opt_algo['Npop'],Nvar))
            # bee(1:OPT_algo.Npop) = struct('pos',zeros(1,Nvar), 'vel',zeros(1,Nvar), 'f',0,'f_p',999999999999999, 'p',zeros(1,Nvar), 'f_g',999999999999999, 'g',zeros(1,Nvar)) 
            for i in range(opt_algo['Npop']):
                bee[i] = {'pos':np.zeros((1,Nvar)), 'vel':np.zeros((1,Nvar)), 'f':0,'f_p':999999999999999, 'p':np.zeros((1,Nvar)), 'f_g':999999999999999, 'g':np.zeros((1,Nvar))}
            
            # Randomize the Values
            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 in range(opt_algo['Npop']):
                for i2 in range(Nvar):
                    bee[i1]['pos'][0,i2] = randomNum(low[i2],high[i2],'dec')
                    bee[i1]['vel'][0,i2] = randomNum(-Vmax[i2],Vmax[i2],'dec')
                    if var['bin'][i2]:
                        bee[i1]['pos'][0,i2] = round(bee[i1]['pos'][0,i2])
            
            # Break the popn into transfer segments
            # ie, popn is 75x40 with 2 transfers, so the first stansfer 
            # segment is popn(:,1:20) and the second segment is popn(:,21:end)
            for i3 in range(opt_algo['Npop']):
                for i4 in range(var['transfers']-1):
                    start = bee[i3]['pos'][0,i4*11]
                    tof = bee[i3]['pos'][0,i4*11+10]
                    bee[i3]['pos'][0,i4*11+11] = start+tof
                
                popn[i3,:] = bee[i3]['pos']
           
        elif algorithm == 'MBH':
            
            # Initials
            bee = 0

            # Pre-Allocate with Random Normalized Values
            popn = np.random.random((opt_algo['N1_Outer'],Nvar)) # Population
            
            # Adjust variables to be within bounds, round binary variables, set
            # constant variables equal to constant
            for i1 in range(opt_algo['N1_Outer']):
                for i2 in range(Nvar):
                    if var['bin'][i2] == 0:      # Non-binary variables
                        popn[i1,i2] = (high[i2]-low[i2])*popn[i1,i2] + low[i2]  
                    elif var['bin'][i2]:   # Binary variables
                        popn[i1,i2] = round(popn[i1,i2])
            
            # Break the popn into transfer segments
            for i3 in range(opt_algo['N1_Outer']):
                for i4 in range(var['transfers']-1):
                    start = popn[i3,i4*11]
                    tof = popn[i3,i4*11+10]
                    popn[i3,i4*11+11] = start+tof
                    
        else:
            
            raise Exception('Invalid Island Algorithm')
        
    else:
        
        raise Exception("Invalid Solver")
    
    return popn,bee