#  FORM: [isBroken,opts] = errorMain(bod,opts)
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function is called from the main script and will abort
# |     execution if any common user errors are input at the start
# |
# |     -The errors are currently:
# |         Variable limits
# |         MBH, percent elements ~= migration elements
# |         Old call to DIR_FSMseg_2D and forcing heliocentric orbits
# |         Optimzation before 1 Jan 1970 or after 31 Dec 2030
# |         
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
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -is_broken          (1,1)       [boolean]   	[unitless]
# |         If any of the user selections will cause the program to break
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import time
import warnings
from julian import to_jd


def errorMain(bod,opts):
    isBroken = False;
    
    
    
# %% Variable Limits 
    
    if len(bod['bodies'])-1 != len(opts['thrust']['tt_end']):
        print('The number of transfers and the list of bodies is not equal.')
        print('Fix "tof_total" or "tof_margin".')
        print('Fix the list of bodies under "bodies".')
        raise Exception()
        
        isBroken = True;
            
    #MBH number of elements in percent selection ~= number of migrations
    
    if opts['MBH']:
                
        if (len(opts['MBH']['per_feas'])) < (opts['island']['Nmig']+1)\
                or (len(opts['MBH']['per_rand'])) < (opts['island']['Nmig']+1):

            print('MBH will break on the last migration.')
            print('The number of migrations is larger than the size of "per_feas" or "per_rand" in "opts[MBH]".\n')
            print('Reduce the number of migrations or add another element to "per_feas" or "per_rand".')
            
            raise Exception()
            isBroken = True;
            return

        elif (len(opts['MBH']['per_feas'])) > (opts['island']['Nmig']+1)\
                or (len(opts['MBH']['per_rand']) > (opts['island']['Nmig']+1)):

            print('MBH has more elements in "per_feas" or "per_rand" than Island Iterations.')
            print("MBH has "+ str(len(opts['MBH']['per_feas'])) + " 'per_feas' parameters and "+str(len(opts['MBH']['per_rand']))+"'per_rand' parameters.")
            print("The optimizer will still run but only the first " + str(opts['island']['Nmig']+1) + " MBH parameter(s) will be used.")
            print('Pausing for 10 sec...')
            time.sleep(10)
                
    # %%Check opts['thrust']
    
    #Having the segmented trajectory have heliocentric orbits within it
    try:
        if opts['thrust']['orbit_check'] =='on':
    
            print('As of right now, MGALT STOpS does not support FORCED multiple heliocentric revolutions.')
            print('This is future work, as the direct solver needs to specify which transfer arc to apply the heliocentric revolution to.')
            
            raise Exception()
            isBroken = True;
    
    except:
        
        #For the indirect methods, which don't have the h-orbit field, make
        #sure to let the warning bypass the error catch
        warnings.warn('Reference to non-existent field')   
    
    
    #If doing MGALT Transfers        
    if opts['solver'] == 'MGALT_DIR_FBSM_2D' or'MGALT_IN_FBSM_2D':
        
       # Check duty cycle
        if bool(opts['thrust'].get('duty_cycle')) == False:
            
            print('For MGALT trajectories, the thruster duty cycle needs to be known.')
            raise Exception()
            isBroken = True;
        
        #Check number of thrusters
        if bool(opts['thrust'].get('n_available')) == False:
            
            print('For MGALT trajectories, the total number of thrusters needs to be known.')
            raise Exception()
            isBroken = True;

    #Prevent the INDIRECT Method from having variable thrust
            
    if opts['solver'] in ['LT_IN_FSM_2D', 'MGALT_IN_FBSM_2D']:
            
        if opts['thrust'].get('thrust_method') == 'variable':
                
            print('The INDIRECT Method does not support varaible thrust.')
            raise Exception()
            isBroken = True;
    
    #Calculate a trajectory before 1 Jan 1970 or past 31 Dec 2030
    
    JD1 = to_jd(bod['window1']) 
    JD2 = to_jd(bod['window2']) + sum(opts['thrust']['tt_end']) + sum(i[1] for i in opts['thrust']['time'])
    
    if JD1 < 2440587.5:
        
        print('MGALT STOpS does not support trajectory optimization before 1 Jan 1970.')
        print('To negate this, new JPL Horizons Data is needed or accurate ode45 support for planet position will have to be implimented.')
        
        raise Exception()
        isBroken = True
    
    if JD2 > 2462867.25:
        
        print('MGALT STOpS does not support trajectory optimization past 31 Dec 2030.')
        print('To negate this, new JPL Horizons Data is needed or accurate ode45 support for planet position will have to be implimented.')
        
        raise Exception()
        isBroken = True
        
    # Make sure Window2/JD2 ~< Window1/JD1
    
    #JD's
    if JD2 < JD1:
        
        print('Julian Date 2 before Julian Date 1, optimization will not work.')
        raise Exception()
        isBroken = True
        
    # Make sure constraint pruning is working
    if opts['constraint_Pruning'] not in ['y','n']:
        raise Exception("Invalid Input for Constraint Pruning On/Off")
    
    # Make sure best leg tracking is working
    if opts['best_Leg_Tracking'] not in ['y','n']:
        raise Exception("Invalid Input for Best Leg Tracking On/Off")
                            
    
    # Launch Windows
    # if JD2 < JD1
        
    #     print('Launch window 2 is before launch window 1, optimization will not work.')
    
    #     isBroken = True
    
    
    
    # User Specified Algorithms ~= Algorithm List

    #Get the island list and the number of those islands from Isl_opt
    GC = [opts['island']['isl_list'].count('MBH'),\
          opts['island']['isl_list'].count('GA'),\
          opts['island']['isl_list'].count('DE'),\
          opts['island']['isl_list'].count('PSO')]
    GR = ['MBH','GA','DE','PSO']
    
    # Run through all of the island lists
    for i1 in range(len(GR)):
        
        # For GA
        if hasattr(GR[i1],'GA'):
            if GC[i1] < len(opts['GA']):
                print('There are more GA islands from "parametersGA" than there are GA connections from "optionsIsland".')
                print("The optimizer will still run but only the first " + GC[i1] + " GA parameter(s) will be used.")
                print('Pausing for 10 sec...')
                time.sleep(10)

            if GC[i1] > len(opts['GA']):
                print('There are more GA connections from "optionsIsland" than there are GA islands from "parametersGA"')
                raise Exception()
                isBroken = True;
    
        
        # For DE
        if hasattr(GR[i1],'DE'):
            if GC[i1] < len(opts['DE']):
                print('There are more DE islands from "parametersDE" than there are DE connections from "optionsIsland".')
                print("The optimizer will still run but only the first " + GC[i1] + "DE parameter(s) will be used.")
                print('Pausing for 10 sec...')
                time.sleep(10)
            
            if GC[i1] > len(opts['DE']):
                print('There are more DE connections from "optionsIsland" than there are DE islands from "parametersDE".')
                raise Exception()
                isBroken = True;
                return
        
        # For PSO
        if hasattr(GR[i1],'PSO'):
            if GC[i1] < len(opts['PSO']):
                print('There are more PSO islands from "parametersPSO" than there are PSO connections from "optionsIsland"')
                print("The optimizer will still run but only the first " + GC[i1] + " PSO parameter(s) will be used.")
                print('Pausing for 10 sec...')
                time.sleep(10)
                
            if GC[i1] > len(opts['PSO']):
                print('There are more PSO connections from "optionsIsland" than there are PSO islands from "parametersPSO".')
                raise Exception()
                isBroken = True;
        
       #  For MBH
        if hasattr(GR[i1],'MBH'):
            if GC[i1] < len(opts['MBH']):
                print('There are more MBH islands from "parametersMBH" than there are MBH connections from "optionsIsland".')
                print("The optimizer will still run but only the first " + GC[i1] + " MBH parameter(s) will be used.")
                print('Pausing for 10 sec...')
                time.sleep(10)
            
            if GC[i1] > len(opts['MBH']):
                print('There are more MBH connections from "optionsIsland" than there are MBH islands from "parametersMBH".')
                raise Exception()
                isBroken = True;    
    
    # Prevent Multiple FSM
    
        
    if opts['solver'] in ["LT_DIR_FSM_2D", "LT_IN_FSM_2D"]:
            
            if len(bod['bodies']) > 2:
                
                print('MGALT STOpS does not support FSM with more than two planets.')
                print('It is recommended to use the FBSM instead.')
                raise Exception()
                isBroken = True;
                return isBroken

    # Prevent Best Leg Tracking on Single Transfer Mission
    
    if len(bod['bodies'])-1 == 1 and opts['best_Leg_Tracking'] == 'y':
        opts['best_Leg_Tracking'] == 'n'
        print('Best Leg Tracking was found to be on. Unfortunately, it cannot be activated for single transfer missions, and is now considered off. Will resume optimization in 10 seconds.')
        time.sleep(10)
        
    if opts['best_Leg_Tracking'] == 'y':
        while opts['island']['isl_list'][0] == 'MBH':
            opts['island']['isl_list'][:-1] = opts['island']['isl_list'][1:]
            opts['island']['isl_list'][-1] = 'MBH'
        

    # Check if Parallel Toolbox
    '''
    try:
    	multiprocessing.Process('distcomp');
        clc
        
    catch
        
        errorPathDisplay();
        fprintf(2,'Parallel Processing Toolbox not installed.\n')
        fprintf(2,'Set "OPT.parallel" to false under the "Main Script" header.\n')
        fprintf(2,'If this message is an error, disable it under "errorMain".\n')
    
        isBroken = true;
        return isBroken
        
    end
    
    
    
    end
    '''
    return isBroken,opts
                
        
        
        
        
        
        
        
