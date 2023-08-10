# PI-MGALT-STOpS Demo - Performance Improved Multiple Gravity Assist Low-Thrust Trajectory Optimization Suite Demo

#%% Good Ol' Imports
from datetime import datetime
from julian import to_jd
import time

from Options.optionsLowThrust import optionsLowThrust
from Options.optionsIsland import optionsIsland
from Options.optionsCost import optionsCost

from Algorithms.Genetic_Algorithm.MGALT_GA import MGALT_GA
from Algorithms.Differential_Evolution.MGALT_DE import MGALT_DE
from Algorithms.Particle_Swarm.MGALT_PSO import MGALT_PSO
from Algorithms.Monotonic_Basin.MGALT_MBH import MGALT_MBH
from Algorithms.Algorithm_Parameters.parametersGA import parametersGA
from Algorithms.Algorithm_Parameters.parametersDE import parametersDE
from Algorithms.Algorithm_Parameters.parametersPSO import parametersPSO
from Algorithms.Algorithm_Parameters.parametersMBH import parametersMBH

from Errors.errorMain import errorMain

from Other_Functions.planetBodies import planetBodies
from Other_Functions.MGALT_varLimits import MGALT_varLimits
from Other_Functions.getConstants import getConstants
from Other_Functions.displayResults import displayResults

#%% Demo Mission
class presets():
    
    def demo_short():
        
        # Initializing...
        bod = {}
        opts = {}
                
        
        #%% Departure/Transfer/Arrival Planet and Earliest/Latest Departure
        """
        Section Description:

        -Each body is a string cooresponding to one of the 9 planets.
        -Options are: 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 
          'Uranus', and 'Neptune'.
        -The position and velocity of the departure body at the start time will 
          be used to create the inital conditions for the trajectory. The position 
          and velocity of the arrival body at the end time will be used to create 
          the end conditions for the trajectory.
        -If the trajectory's goal is not to arrive at the arrival body, but just 
          to insert into its orbit, the arrival body's location will be ignored 
          and the end conditions become the position and velocity at some point 
          on its orbit.
        -All entries are spelling and case sensitive
        """
    
        bod['bodies'] = ('Earth','Mars','Jupiter')
        
        #%% Launch Window
        
        """
        Section Description:
            
        -Here the user is required to enter the earliest and latest departure date
          from the departure body. The earliest date is refered to as window1 and
          the latest departure date is refered to as window2.  The format for the
          date is datetime(YEAR,MONTH,DAY) with each input (year, month, and day) being an
          integer.
        """
        bod['window1'] = datetime(2020,11,1)       # earliest launch day	[year,month,day]
        bod['window2'] = datetime(2022,12,1)       # latest launch day  	[year,month,day]
        
        #%% Time of Flight Options
        """
        Section Description:
            
        -Time of flight parameters for the whole solver in days.
          It is important to specify the desired transfer time between each
          transfer segment, with the allowable tolerance.
            -Ex: On a mission from Earth->Mars->Jupiter, there are 2 flight
            segments which need a time specified for them
                tof_total = [(300),(900)]
                tof_margin = [(50,50)],(250,250)]
            This tells the solver that a trajectory time of 300-+50 days is desired
            for the transfer from Earth->Mars, and 900-+250 days is desired from
            Mars->Jupiter
        -Add margin to the desired days for random perturbations
        """
        opts['tof_total'] = [(752),\
                             (1448)] # tof w/ gravity assists [days]
            
        opts['tof_margin'] = [(100,100),\
                              (300,300)] # margin on tof -/+     [days]    
    
        #%% Solver and Thrust Method
        """
        Section Description:

        -The chosen cost funtion dictates how the trajectory is described. A 
          variable string decribes each possible trajectory. In this work there 
          are two main structures for the variable strings.

        -In what will be called the segmented method the trajectory is divided 
          into N segments each with a thrust value and a thrust pointing angle. 
          The variable string for the segmented method consists of a departure 
          time, the thrust and pointing angle for each segment, and an arrival 
          time. 

        -In what will be called the costate method the trajectory has a consant 
          or equation thrust and the time dependent thrust angle is described 
          with three intial costate variables. The variable string for the 
          Conway method consists of a departure time, the three costate variables, 
          and an arrival time. 

        -The two cost function handles available coorespond to the two variable 
          composition methods. They are: 'EP_cost_fun_segmented_2D' and 
          'EP_cost_fun_costate_2D'. The cost function is spelling and case 
          sensitive.

        """
        opts['solver'] = 'MGALT_IN_FBSM_2D'
        
        #%% Optimization Options
        
        opts['thrust'] = optionsLowThrust(opts,'Yam-STOUR')
        
            #%%% Island Model
        
        opts['island'] = optionsIsland('2M-all')
        
            #%%% Cost Function
        
        [opts['cost'],opts['weighting'],opts['ode']] = optionsCost(opts,'default')

            #%%% Algorithm Parameters
        
        opts['GA'] = parametersGA('75_30')#     # Genetic Algorithm Island #1
        opts['DE'] = parametersDE('75_30')#      # Differential Evolution Island #1
        opts['PSO'] = parametersPSO('50_50')#     # Particle Swarm Island #1
        opts['MBH'] = parametersMBH('1450_75')#     # Particle Swarm Island

            #%%% Additional Features 
        opts['parallel'] = 'y' # 'y' or 'n' #Decides Parallel Programming
        opts['constraint_Pruning'] = 'n' # 'y' or 'n' # Recommended on for Direct Method, off for Indirect Method
        opts['cluster_Pruning'] = 'y' # 'y' or 'n' # Recommended on for longer runs, off for short runs
        opts['best_Leg_Tracking'] = 'y'  # 'y' or 'n' # Decides Best Leg Solution Tracking
        
        return bod, opts
        
    def demo_long():
        
        # Initializing...
        bod = {}
        opts = {}
                
        
        #%% Departure/Transfer/Arrival Planet and Earliest/Latest Departure
        """
        Section Description:

        -Each body is a string cooresponding to one of the 9 planets.
        -Options are: 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 
          'Uranus', and 'Neptune'.
        -The position and velocity of the departure body at the start time will 
          be used to create the inital conditions for the trajectory. The position 
          and velocity of the arrival body at the end time will be used to create 
          the end conditions for the trajectory.
        -If the trajectory's goal is not to arrive at the arrival body, but just 
          to insert into its orbit, the arrival body's location will be ignored 
          and the end conditions become the position and velocity at some point 
          on its orbit.
        -All entries are spelling and case sensitive
        """
    
        bod['bodies'] = ('Earth','Mars','Jupiter')
        
        #%% Launch Window
        
        """
        Section Description:
            
        -Here the user is required to enter the earliest and latest departure date
          from the departure body. The earliest date is refered to as window1 and
          the latest departure date is refered to as window2.  The format for the
          date is datetime(YEAR,MONTH,DAY) with each input (year, month, and day) being an
          integer.
        """
        bod['window1'] = datetime(2020,11,1)       # earliest launch day	[year,month,day]
        bod['window2'] = datetime(2022,12,1)       # latest launch day  	[year,month,day]
        
        #%% Time of Flight Options
        """
        Section Description:
            
        -Time of flight parameters for the whole solver in days.
          It is important to specify the desired transfer time between each
          transfer segment, with the allowable tolerance.
            -Ex: On a mission from Earth->Mars->Jupiter, there are 2 flight
            segments which need a time specified for them
                tof_total = [(300),(900)]
                tof_margin = [(50,50)],(250,250)]
            This tells the solver that a trajectory time of 300-+50 days is desired
            for the transfer from Earth->Mars, and 900-+250 days is desired from
            Mars->Jupiter
        -Add margin to the desired days for random perturbations
        """
        opts['tof_total'] = [(752),\
                             (1448)] # tof w/ gravity assists [days]
            
        opts['tof_margin'] = [(100,100),\
                              (300,300)] # margin on tof -/+     [days]    
    
        #%% Solver and Thrust Method
        """
        Section Description:

        -The chosen cost funtion dictates how the trajectory is described. A 
          variable string decribes each possible trajectory. In this work there 
          are two main structures for the variable strings.

        -In what will be called the segmented method the trajectory is divided 
          into N segments each with a thrust value and a thrust pointing angle. 
          The variable string for the segmented method consists of a departure 
          time, the thrust and pointing angle for each segment, and an arrival 
          time. 

        -In what will be called the costate method the trajectory has a consant 
          or equation thrust and the time dependent thrust angle is described 
          with three intial costate variables. The variable string for the 
          Conway method consists of a departure time, the three costate variables, 
          and an arrival time. 

        -The two cost function handles available coorespond to the two variable 
          composition methods. They are: 'EP_cost_fun_segmented_2D' and 
          'EP_cost_fun_costate_2D'. The cost function is spelling and case 
          sensitive.

        """
        opts['solver'] = 'MGALT_IN_FBSM_2D'
        
        #%% Optimization Options
        
        opts['thrust'] = optionsLowThrust(opts,'Yam-STOUR')
        
            #%%% Island Model
        
        opts['island'] = optionsIsland('4M-all')
        
            #%%% Cost Function
        
        [opts['cost'],opts['weighting'],opts['ode']] = optionsCost(opts,'default')

            #%%% Algorithm Parameters
        
        opts['GA'] = parametersGA('450_180')#     # Genetic Algorithm Island
        opts['DE'] = parametersDE('450_180')#      # Differential Evolution Island
        opts['PSO'] = parametersPSO('350_150')#     # Particle Swarm Island
        opts['MBH'] = parametersMBH('15000_1500')#     # Particle Swarm Island

            #%%% Additional Features 
        opts['parallel'] = 'y' # 'y' or 'n' #Decides Parallel Programming
        opts['constraint_Pruning'] = 'n' # 'y' or 'n' # Recommended on for Direct Method, off for Indirect Method
        opts['cluster_Pruning'] = 'y' # 'y' or 'n' # Recommended on for longer runs, off for short runs
        opts['best_Leg_Tracking'] = 'y'  # 'y' or 'n' # Decides Best Leg Solution Tracking
        
        return bod, opts
    
    # Mission Initialization Complete...
    
#%% Run Mission
def run(miss):
    # # # Core Functionality # # #
    start = time.time() #start timer

    if miss in dir(presets):
        func = getattr(presets, miss)
        [bod,opts] = func()
    else:
        raise Exception("Please enter valid mission name.")
    
    [isBroken,opts] = errorMain(bod,opts)
    
    opts['ephemeris'] = 'JPL_Horizons' #or 'jplEphem'
    bod['windows'] = [bod['window1'],bod['window2']]
    bod['JD1'] = to_jd(bod['windows'][0])
    bod['JD2'] = to_jd(bod['windows'][1]) + sum(opts['thrust']['tt_end']) + sum(i[1] for i in opts['thrust']['time'])
    bod['tspan'] = [bod['JD1'],bod['JD2']]
    bod['bodies_R'],bod['bodies_V'],bod['bodies_JD'] = planetBodies(bod,opts)
    var = MGALT_varLimits(bod,opts)
    const = getConstants(bod)
    
    #%% Optimization

        # Section Description:
    '''
    -Here is where the optimization is actually run. All of the user inputs
      from the previous sections are run through here into the optimization
      algorithms.
    '''
    selected = [[{'f':[],'sol':[]} for iter in range(len(opts['island']['isl_list']))] for it in range(opts['island']['Nmig']+1)]
    eval_info = [[{'optimal_soln':[],'f_best':[],'iterations':[],'maxcost':[],'mincost':[],'avgcost':[],'total_evals':[]} for iter in range(len(opts['island']['isl_list']))] for it in range(opts['island']['Nmig']+1)]

    # Run through the solver
    for num_mig in range(opts['island']['Nmig']+1):

        # Algorithm counter in case an algorithm is used more than once
        count_GA = 0 
        count_DE = 0 
        count_PSO = 0 
        count_MBH = 0
        
        # Run through all the islands
        for num_isl in range(opts['island']['Nisl']):

            if opts['island']['isl_list'][num_isl] == 'GA':
            
                [eval_info[num_mig][num_isl], selected,fs,fs_trans_GA] = MGALT_GA(bod,const,opts,var,selected,num_mig,num_isl,count_GA)
                count_GA = count_GA + 1
                
                bees = 0
            
            elif opts['island']['isl_list'][num_isl] == 'DE':
                
                [eval_info[num_mig][num_isl], selected,fs,fs_trans_DE] = MGALT_DE(bod,const,opts,var,selected,num_mig,num_isl,count_DE) #ok<*SAGROW>
                count_DE = count_DE + 1

                bees = 0
                
            elif opts['island']['isl_list'][num_isl] == 'PSO':
                
                [eval_info[num_mig][num_isl], selected,fs,bees,fs_trans_PSO] = MGALT_PSO(bod,const,opts,var,selected,num_mig,num_isl,count_PSO) #ok<*SAGROW>
                count_PSO = count_PSO + 1
                
            elif opts['island']['isl_list'][num_isl] == 'MBH':
                
                [eval_info[num_mig][num_isl], selected,fs,fs_trans_MBH] = MGALT_MBH(bod,const,opts,var,selected,num_mig,num_isl,count_MBH) #ok<*SAGROW>
                count_MBH = count_MBH + 1
                
                bees = 0
                
            else:
                
                raise Exception("Invalid Island")
                
        #%% Gather Results
    print('Migrations Complete!')
    end = time.time() #end timer
    run_time = end-start
    best_sol = displayResults(bod,const,opts,var,run_time,eval_info)
    return bod,opts,var,isBroken,best_sol