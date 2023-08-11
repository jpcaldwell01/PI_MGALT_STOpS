# PI-MGALT-STOpS - Performance Improved Multiple Gravity Assist Low-Thrust Trajectory Optimization Suite Main Script

#%% Good Ol'Imports
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
from Other_Functions.plotPareto import plotPareto

#%% Presets

# For instructions on creating a preset, see STOpS_Demo.py

class presets():
    
    def custom():
        
        bod = {}
        opts = {}
    
        bod['bodies'] = ('Earth','Mars')
        opts['tof_total'] = [(752)] # tof w/ gravity assists [days]
        opts['tof_margin'] = [(100,100)] # margin on tof -/+     [days]

        # Launch Windows
        bod['window1'] = datetime(2021,10,1)       # earliest launch day	[year,month,day]
        bod['window2'] = datetime(2022,12,1)       # latest launch day  	[year,month,day]
        
        # Solver
        # opts['solver'] = 'LT_DIR_FSM_2D'
        # opts['solver'] = 'LT_IN_FSM_2D'
        # opts['solver'] = 'MGALT_DIR_FBSM_2D'
        opts['solver'] = 'MGALT_IN_FBSM_2D'
        
        opts['thrust'] = optionsLowThrust(opts,'Yam-STOUR')
        opts['island'] = optionsIsland('2M-all')
        [opts['cost'],opts['weighting'],opts['ode']] = optionsCost(opts,'default')
        
        # Algorithm Parameters
        opts['GA'] = parametersGA('75_30')#     # Genetic Algorithm Island #1
        opts['DE'] = parametersDE('75_30')#      # Differential Evolution Island #1
        opts['PSO'] = parametersPSO('50_50')#     # Particle Swarm Island #1
        opts['MBH'] = parametersMBH('1450_75')#     # Particle Swarm Island
        
        # Performance Improvements
        opts['parallel'] = 'n' # 'y' or 'n' #Decides Parallel Programming
        opts['constraint_Pruning'] = 'n' # 'y' or 'n' # Recommended on for Direct Method, off for Indirect Method
        opts['cluster_Pruning'] = 'n' # 'y' or 'n' # Recommended on for longer runs, off for short runs
        opts['best_Leg_Tracking'] = 'y'  # 'y' or 'n' # Decides Best Leg Solution Tracking
        
        return bod, opts
    
    def EVMVEJ_IN():
        
        bod = {}
        opts = {}
    
        bod['bodies'] = ('Earth','Venus','Mars','Venus','Earth','Jupiter')
        opts['tof_total'] = [(150),\
                             (128),\
                             (330),\
                             (58),\
                             (770)] # tof w/ gravity assists [days]
        opts['tof_margin'] = [(50,50),\
                              (50,50),\
                              (50,50),\
                              (15,50),\
                              (50,50)] # margin on tof -/+     [days]

        # Launch Windows
        bod['window1'] = datetime(2021,9,24)       # earliest launch day	[year,month,day]
        bod['window2'] = datetime(2021,11,24)       # latest launch day  	[year,month,day]
        
        # Solver
        # opts['solver'] = 'LT_DIR_FSM_2D'
        # opts['solver'] = 'LT_IN_FSM_2D'
        # opts['solver'] = 'MGALT_DIR_FBSM_2D'
        opts['solver'] = 'MGALT_IN_FBSM_2D'
        
        opts['thrust'] = optionsLowThrust(opts,'Yam-STOUR')
        opts['island'] = optionsIsland('4M-all')
        [opts['cost'],opts['weighting'],opts['ode']] = optionsCost(opts,'default')
        
        # Algorithm Parameters
        opts['GA'] = parametersGA('450_180')#     # Genetic Algorithm Island #1
        opts['DE'] = parametersDE('450_180')#      # Differential Evolution Island #1
        opts['PSO'] = parametersPSO('350_150')#     # Particle Swarm Island #1
        opts['MBH'] = parametersMBH('15000_1500')#     # Particle Swarm Island
        
        # Performance Improvements
        opts['parallel'] = 'y' # 'y' or 'n' # Decides Parallel Programming
        opts['constraint_Pruning'] = 'n' # 'y' or 'n' # Recommended on for Direct Method, off for Indirect Method
        opts['cluster_Pruning'] = 'y' # 'y' or 'n' # Recommended on for longer runs, off for short runs
        opts['best_Leg_Tracking'] = 'y'  # 'y' or 'n' # Decides Best Leg Solution Tracking
        
        return bod, opts
    
    def MGALT_STOpS_T2_IN_SHORT():
        
        bod = {}
        opts = {}
    
        bod['bodies'] = ('Earth','Mars','Jupiter')
        opts['tof_total'] = [(752),\
                             (1448)] # tof w/ gravity assists [days]
        opts['tof_margin'] = [(100,100),\
                              (300,300)] # margin on tof -/+     [days]

        # Launch Windows
        bod['window1'] = datetime(2020,11,1)       # earliest launch day	[year,month,day]
        bod['window2'] = datetime(2022,12,1)       # latest launch day  	[year,month,day]
        
        # Solver
        # opts['solver'] = 'LT_DIR_FSM_2D'
        # opts['solver'] = 'LT_IN_FSM_2D'
        # opts['solver'] = 'MGALT_DIR_FBSM_2D'
        opts['solver'] = 'MGALT_IN_FBSM_2D'
        
        opts['thrust'] = optionsLowThrust(opts,'Yam-STOUR')
        opts['island'] = optionsIsland('2M-all')
        [opts['cost'],opts['weighting'],opts['ode']] = optionsCost(opts,'default')
        
        # Algorithm Parameters
        opts['GA'] = parametersGA('75_30')#     # Genetic Algorithm Island #1
        opts['DE'] = parametersDE('75_30')#      # Differential Evolution Island #1
        opts['PSO'] = parametersPSO('50_50')#     # Particle Swarm Island #1
        opts['MBH'] = parametersMBH('1450_75')#     # Particle Swarm Island
        
        # Performance Improvements
        opts['parallel'] = 'n' # 'y' or 'n' #Decides Parallel Programming
        opts['constraint_Pruning'] = 'y' # 'y' or 'n' # Recommended on for Direct Method, must turn off for Indirect Method
        opts['cluster_Pruning'] = 'y' # 'y' or 'n' # Recommended on for longer runs, off for short runs
        opts['best_Leg_Tracking'] = 'y'  # 'y' or 'n' # Decides Best Leg Solution Tracking


        return bod, opts
    
    def MGALT_STOpS_T2_IN_LONG():
        
        bod = {}
        opts = {}
    
        bod['bodies'] = ('Earth','Mars','Jupiter')
        opts['tof_total'] = [(752),\
                             (1448)] # tof w/ gravity assists [days]
        opts['tof_margin'] = [(100,100),\
                              (300,300)] # margin on tof -/+     [days]

        # Launch Windows
        bod['window1'] = datetime(2020,11,1)       # earliest launch day	[year,month,day]
        bod['window2'] = datetime(2022,12,1)       # latest launch day  	[year,month,day]
        
        # Solver
        # opts['solver'] = 'LT_DIR_FSM_2D'
        # opts['solver'] = 'LT_IN_FSM_2D'
        # opts['solver'] = 'MGALT_DIR_FBSM_2D'
        opts['solver'] = 'MGALT_IN_FBSM_2D'
        
        opts['thrust'] = optionsLowThrust(opts,'Yam-STOUR')
        opts['island'] = optionsIsland('4M-all')
        [opts['cost'],opts['weighting'],opts['ode']] = optionsCost(opts,'default')
        
        # Algorithm Parameters
        opts['GA'] = parametersGA('450_180')#     # Genetic Algorithm Island #1
        opts['DE'] = parametersDE('450_180')#      # Differential Evolution Island #1
        opts['PSO'] = parametersPSO('350_150')#     # Particle Swarm Island #1
        opts['MBH'] = parametersMBH('15000_1500')#     # Particle Swarm Island
        
        # Performance Improvements
        opts['parallel'] = 'y' # 'y' or 'n' # Decides Parallel Programming
        opts['constraint_Pruning'] = 'n' # 'y' or 'n' # Recommended on for Direct Method, off for Indirect Method
        opts['cluster_Pruning'] = 'y' # 'y' or 'n' # Recommended on for longer runs, off for short runs
        opts['best_Leg_Tracking'] = 'y'  # 'y' or 'n' # Decides Best Leg Solution Tracking


        return bod, opts

    def MGALT_STOpS_T2_DIR_SHORT():
        
        bod = {}
        opts = {}
    
        bod['bodies'] = ('Earth','Mars','Jupiter')
        opts['tof_total'] = [(855),\
                             (1346)] # tof w/ gravity assists [days]
        opts['tof_margin'] = [(50,50),\
                              (50,50)] # margin on tof -/+     [days]

        # Launch Windows
        bod['window1'] = datetime(2021,9,16)       # earliest launch day	[year,month,day]
        bod['window2'] = datetime(2022,1,16)       # latest launch day  	[year,month,day]
        
        # Solver
        # opts['solver'] = 'LT_DIR_FSM_2D'
        # opts['solver'] = 'LT_IN_FSM_2D'
        opts['solver'] = 'MGALT_DIR_FBSM_2D'
        # opts['solver'] = 'MGALT_IN_FBSM_2D'
        
        opts['thrust'] = optionsLowThrust(opts,'Yam-GALLOP')
        opts['island'] = optionsIsland('2M-all')
        [opts['cost'],opts['weighting'],opts['ode']] = optionsCost(opts,'default')
        
        # Algorithm Parameters
        opts['GA'] = parametersGA('75_30')#     # Genetic Algorithm Island #1
        opts['DE'] = parametersDE('75_30')#      # Differential Evolution Island #1
        opts['PSO'] = parametersPSO('50_50')#     # Particle Swarm Island #1
        opts['MBH'] = parametersMBH('750_300')#     # Particle Swarm Island
        
        # Performance Improvements
        opts['parallel'] = 'y' # 'y' or 'n' #Decides Parallel Programming
        opts['constraint_Pruning'] = 'y' # 'y' or 'n' # Recommended on for Direct Method, off for Indirect Method
        opts['cluster_Pruning'] = 'y' # 'y' or 'n' # Recommended on for longer runs, off for short runs
        opts['best_Leg_Tracking'] = 'y'  # 'y' or 'n' # Decides Best Leg Solution Tracking


        return bod, opts
    
    def MGALT_STOpS_T2_DIR_LONG():
        
        bod = {}
        opts = {}
    
        bod['bodies'] = ('Earth','Mars','Jupiter')
        opts['tof_total'] = [(855),\
                             (1346)] # tof w/ gravity assists [days]
        opts['tof_margin'] = [(50,50),\
                              (50,50)] # margin on tof -/+     [days]

        # Launch Windows
        bod['window1'] = datetime(2021,8,1)       # earliest launch day	[year,month,day]
        bod['window2'] = datetime(2022,3,1)       # latest launch day  	[year,month,day]
         
        # Solver
        # opts['solver'] = 'LT_DIR_FSM_2D'
        # opts['solver'] = 'LT_IN_FSM_2D'
        opts['solver'] = 'MGALT_DIR_FBSM_2D'
        # opts['solver'] = 'MGALT_IN_FBSM_2D'
        
        opts['thrust'] = optionsLowThrust(opts,'Yam-GALLOP')
        opts['island'] = optionsIsland('3M-all')
        [opts['cost'],opts['weighting'],opts['ode']] = optionsCost(opts,'default')
        
        # Algorithm Parameters
        opts['GA'] = parametersGA('300_115')#     # Genetic Algorithm Island
        opts['DE'] = parametersDE('300_115')#      # Differential Evolution Island
        opts['PSO'] = parametersPSO('250_100')#     # Particle Swarm Island
        opts['MBH'] = parametersMBH('10000_800')#     # Particle Swarm Island
        
        # Performance Improvements
        opts['parallel'] = 'y'           # 'y' or 'n' # Decides Parallel Programming
        opts['constraint_Pruning'] = 'y' # 'y' or 'n' # Recommended on for Direct Method, off for Indirect Method
        opts['cluster_Pruning'] = 'y'    # 'y' or 'n' # Recommended on for longer runs, off for short runs
        opts['best_Leg_Tracking'] = 'y'  # 'y' or 'n' # Decides Best Leg Solution Tracking


        return bod, opts

#%% Run Mission

def run(miss):
    
    # # # Core Functionality # # #
    start = time.time() #start timer

    if miss in dir(presets):
        func = getattr(presets, miss)
        [bod,opts] = func()
    else:
        raise Exception("Please enter valid mission name.")
    
    opts['ephemeris'] = 'JPL_Horizons'
    [isBroken,opts] = errorMain(bod,opts)
    
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
            
            elif opts['island']['isl_list'][num_isl] == 'DE':
                
                [eval_info[num_mig][num_isl], selected,fs,fs_trans_DE] = MGALT_DE(bod,const,opts,var,selected,num_mig,num_isl,count_DE) #ok<*SAGROW>
                count_DE = count_DE + 1
                
            elif opts['island']['isl_list'][num_isl] == 'PSO':
                
                [eval_info[num_mig][num_isl], selected,fs,bees,fs_trans_PSO] = MGALT_PSO(bod,const,opts,var,selected,num_mig,num_isl,count_PSO) #ok<*SAGROW>
                count_PSO = count_PSO + 1
                
            elif opts['island']['isl_list'][num_isl] == 'MBH':
                
                [eval_info[num_mig][num_isl], selected,fs,fs_trans_MBH] = MGALT_MBH(bod,const,opts,var,selected,num_mig,num_isl,count_MBH) #ok<*SAGROW>
                count_MBH = count_MBH + 1
                
            else:
                
                raise Exception("Invalid Island")
                
        #Plotting Migration Pareto
        plot_pareto = 1 # Change to 1 to see 2D Pareto Plot. Built for 2 transfer trajectory
        if num_mig == opts['island']['Nmig'] and len(bod['bodies']) > 2 and plot_pareto == 1:
            plotPareto(fs_trans_GA,fs_trans_DE,fs_trans_PSO,fs_trans_MBH,num_mig,opts)
                
    #%% Gather Results
        
    print('Migrations Complete!')
    end = time.time() #end timer
    run_time = end-start
    best_sol = displayResults(bod,const,opts,var,run_time,eval_info)
    return bod,opts,var,isBroken,best_sol