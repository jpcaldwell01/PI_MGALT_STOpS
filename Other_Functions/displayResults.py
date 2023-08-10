# FORM: [] = displayResults(bod,const,opts,var,run_time,eval_info)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Function designed to display optimization results
# |
# |     -This funciton gives info regarding optimization time and 
# |     iterations, orbital transfers, and other information regarding 
# |     some of the inputs used for STOpS
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
# |     -VAR                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |     -run_time        	(1,1)       [float]        	[min]
# |         How long the optimization process took to run
# |     -eval_info          (1,1)       [struct]        [unitless]
# |         Best solutions, best costs, number of iterations, etc...
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------
#%%

import numpy as np
from Solvers.Indirect import MGALT_IN_FBSM_2D
from Solvers.Direct import MGALT_DIR_FBSM_2D
from IPython import get_ipython
from julian import from_jd
from Other_Functions.plotOrbits import plotOrbits
#%%
def displayResults(bod,const,opts,var,run_time,eval_info):
    
    ## Plot Results
    
    #Section Description:
    '''
    -Check through all of the results in order to fund the best results for the
      optimization methods.
    '''
    
    # ***Calculated***
    final_soln = np.zeros((opts['island']['Nisl'],len(var['bin'])))
    final_f = np.zeros((1,opts['island']['Nisl']))
    
    for i3 in range(opts['island']['Nisl']):
        final_soln[i3,:] = eval_info[opts['island']['Nmig']][i3]['optimal_soln'][0,:]
        final_f[0,i3] = eval_info[opts['island']['Nmig']][i3]['f_best'][0][0]
    
    I = np.argsort(final_f)
    sorted_final_soln = final_soln[I][0,:,:]
    
    # Plot the results
    if opts['solver'] ==  'MGALT_IN_FBSM_2D':
        feval = getattr(MGALT_IN_FBSM_2D,'MGALT_IN_FBSM_2D')
    elif opts['solver'] ==  'MGALT_DIR_FBSM_2D':
        feval = getattr(MGALT_DIR_FBSM_2D,'MGALT_DIR_FBSM_2D')

    [f,plot_vars,f_trans] = feval(sorted_final_soln[0,:],bod,const,opts,var)
    plotOrbits(bod,const,opts,var,plot_vars)
    run_time = run_time/60
    run_time_hrs = run_time/60.
    ## Display Completion Message
    
    get_ipython().magic('clear')
    print('---------------------\n')
    print('Optimization Complete\n')
    print('---------------------\n\n')
    print('   Solver Used: ',opts['solver'],'\n\n')
    print('Total Run Time: %.2f minutes' % run_time)
    print('                %.2f hours\n' % run_time_hrs)
    print('      Parallel: ',opts['parallel'],' \n')
    print('      Constraint Pruning: ',opts['constraint_Pruning'],' \n')
    print('      Cluster Pruning: ',opts['cluster_Pruning'],' \n')
    print('      Leg Tracking: ',opts['best_Leg_Tracking'],' \n')

    if opts['cluster_Pruning'] == 'y':
        print('          Cluster Pruning Results:')
        if len(var['clusterResults']) > 0:
            date1 = from_jd(var['clusterResults'][0,0])
            date2 = from_jd(var['clusterResults'][0,1])

            print('               Departure Window from',bod['bodies'][0],': ',date1.strftime("%d/%m/%Y"),'-->',date2.strftime("%d/%m/%Y"),'   (DD/MM/YYYY)')
            
            for it in range(var['transfers']):
                print('               ToF from',bod['bodies'][it], 'to',bod['bodies'][it+1],': %.2f --> %.2f days' % (var['clusterResults'][it+1,0],var['clusterResults'][it+1,1]))
            print('\n')
        else:
            print('               No significant clusters found. It is recommended to try larger departure and ToF windows, as well as change the parameters of the DE algorithm to include more evaluations.')
    
    total = 0
    for i1 in range(np.shape(eval_info)[0]):
        for i2 in range(np.shape(eval_info)[1]):
            total = total+eval_info[i1][i2]['total_evals']
    
    print('   Total Iters: ',total,'\n')
    print('     Best Total Cost: %.2f \n' % f)
    if var['transfers'] > 1:
        for i in range(var['transfers']):
            print('        Cost from',bod['bodies'][i],'to',bod['bodies'][i+1],': %f' % f_trans[0,i])
    print('\n\n\n')
    
    ## Orbit Information
    
    TU = const['TU'][0]
    
    print('-----------------\n')
    print('Orbit Information\n')
    print('-----------------\n\n')
    print('                *  DD/MM/YYYY  HH:MM:SS  *\n\n')
    
    date_info = from_jd(sorted_final_soln[0,0])
    print('   Departure Body: ',bod['bodies'][0],'\n')
    print('   Departure Date: ',date_info.strftime("%d/%m/%Y"),' ',date_info.strftime("%H:%M:%S"))
    print('\n====>\n\n')
    
    ToF = np.zeros((var['transfers']-1,1))
    # Display the transfers
    if var['transfers'] == 1:
        
        # Time fo Flight
        if opts['solver'] == 'LT_DIR_FSM_2D':
            1/0
            # ToF = plot_vars.tspan{1,1}(end,end)*TU
        elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
            ToF = plot_vars['tspan_fs'][-1,-1]*2*TU
        elif opts['solver'] in ['LT_IN_FSM_2D','MGALT_IN_FBSM_2D']:
            ToF = plot_vars['tspan'][-1,-1]*TU

        # Arrival Date
        date = from_jd(plot_vars['JD'][-1,-1])

        # Fuel Usage
        data = plot_vars['transfers'][:,-1]
        fuel_start = data[0]
        fuel_end = data[-1]
        
        print('     Arrival Body: ',bod['bodies'][-1],'\n')
        print('   Time of Flight: %.2f days\n' % ToF)
        print('     Arrival Date: ',date_info.strftime("%d/%m/%Y"),' ',date_info.strftime("%H:%M:%S"))
        print('   Fuel Mass Start: %.2f kg\n' % fuel_start)
        print('     Fuel Mass End: %.2f kg\n' % fuel_end)

        # Fuel warning
        if fuel_end <= 1:
            print('\nWARNING!\n')
            print('Fuel dropped below 1kg. \n')
            print('Most likely the orbit did not converge with',bod['bodies'][-1],'\n')
            print('Re-run the solver with different engine parameters!\n')

        print('\n. . . . . . . . . . . . . . . .\n\n')
        
    else:
        
        for i4 in range(1,var['transfers']):

            # Time fo Flight
            if opts['solver'] == 'LT_DIR_FSM_2D':
                    ToF[i4-1] = plot_vars['tspan'][i4-1,-1][-1,-1]*TU
            elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
                    ToF[i4-1] = plot_vars['tspan_fs'][i4-1][-1,-1]*2*TU
            if opts['solver'] in ['LT_IN_FSM_2D','MGALT_IN_FBSM_2D']:
                    ToF[i4-1] = plot_vars['tspan'][i4-1,-1]*TU

            # Flyby Dates
            flyby = plot_vars['JD'][i4-1,-1]
            date_info = from_jd(flyby)

            # Fuel Usage
            data = plot_vars['transfers'][:,(i4-1)*int(np.shape(plot_vars['transfers'])[1]/var['transfers']):i4*int(np.shape(plot_vars['transfers'])[1]/var['transfers'])]
            fuel_start = data[0,-1]
            fuel_end = data[-1,-1]

            print('  Transfer Body',i4,':',bod['bodies'][i4],'\n')
            print('    Time of Flight: %.2f days\n' % ToF[i4-1])
            print('        Flyby Date:',date_info.strftime("%d/%m/%Y"),'',date_info.strftime("%H:%M:%S"),'\n')
            print('   Fuel Mass Start: %.2f kg\n' % fuel_start)
            print('     Fuel Mass End: %.2f kg\n' % fuel_end)

            # Fuel warning
            if fuel_end <= 0.5:
                print('\nWARNING!\n')
                print('Fuel dropped below 0.5kg on Transfer Segment',i4,'\n')
                print('Most likely the orbit did not converge with',bod['bodies'][i4],'\n')
                print('Re-run the solver with different engine parameters!\n\n')
                print('This is just a warning and not an actual error, nothing was suspended during operation.\n')

            print('\n====>\n\n')

        # Time of Flight
        if opts['solver'] == 'LT_DIR_FSM_2D':
                ToF = np.append(ToF,plot_vars['tspan'][-1,-1][-1,-1]*TU)
        elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
                ToF = np.append(ToF,plot_vars['tspan_fs'][-1][-1,-1]*2*TU)
        elif opts['solver'] in ['LT_IN_FSM_2D','MGALT_IN_FBSM_2D']:
                ToF = np.append(ToF,plot_vars['tspan'][-1,-1]*TU)

        # Arrival Date
        date_info = from_jd(plot_vars['JD'][-1,-1])

        # Fuel Usage
        data = plot_vars['transfers'][:,int(np.shape(plot_vars['transfers'])[1]*(var['transfers']-1)/var['transfers']):]
        fuel_start = data[0,-1]
        fuel_end = data[-1,-1]

        print('     Arrival Body: ',bod['bodies'][-1],'\n')
        print('   Time of Flight: %.2f days\n' % ToF[-1])
        print('     Arrival Date: ',date_info.strftime("%d/%m/%Y"),' ',date_info.strftime("%H:%M:%S"))
        print('   Fuel Mass Start: %.2f kg\n' % fuel_start)
        print('     Fuel Mass End: %.2f kg\n' % fuel_end)

        # Fuel warning
        if fuel_end <= 1:
            print('\nWARNING!\n')
            print('Fuel dropped below 1kg on Transfer Segment',i4+1,'\n')
            print('Most likely the orbit did not converge with',bod['bodies'][-1],'\n')
            print('Re-run the solver with different engine parameters!\n')

        print('\n. . . . . . . . . . . . . . . .\n\n')
    
    # Print off the total flight time
    if var['transfers'] == 1:
        print('Total Flight Time: %.2f days\n' % ToF)
        print('                 : %.2f years\n\n' % (ToF/365.25))
    else:
        print('Total Flight Time: %.2f days\n' % np.sum(ToF,axis=None))
        print('                 : %.2f years\n\n' % (np.sum(ToF)/365.25))
    
    #Print off the fuel mass fraction
    mass_total = opts['thrust']['m0']
    mass_structure = plot_vars['transfers'][-1,-1]
    mass_fuel = mass_total-mass_structure
    mass_fraction = mass_fuel/mass_total * 100
    print('   Prop mass frac: %.2f percent\n' % mass_fraction)
    
    print('\n\n\n')
    
    
    
    ## Other Information
    
    print('-----------------\n')
    print('Other Information\n')
    print('-----------------\n\n')
    
    # Search Window
    print('    Search Window:\n')
    print('    ',bod['window1'].strftime("%Y %m %d"),'\n')
    print('    ',bod['window2'].strftime("%Y %m %d"),'\n')
    print('\n. . . . . . . . . . . . . . . .\n\n')
    
    # ToF Parameters
    print('    Desired ToF:\n\n')
    for i5 in range(var['transfers']):
       print('    ',bod['bodies'][i5],' to ',bod['bodies'][i5+1],'\n')
       print('    %.0f Days\n' % opts['thrust']['tt_end'][i5])
       print('    -%4.0f +%4.0f Days\n' % (opts['thrust']['time'][i5][0],opts['thrust']['time'][i5][1]))
       print('\n')

    print('\n. . . . . . . . . . . . . . . .\n\n')
    
    # Thrust Parameters
    print('    Engine Parameters:\n')
    print('    ',opts['thrust']['thrust_method'],'\n')
    if opts['solver'] == 'MGALT_IN_FBSM_2D': 
        print('    @ %2.4f N\n' % opts['thrust']['thrust'])
    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
        print('    @ %2.4f - %2.4f N\n' % (opts['thrust']['thrust'][0],opts['thrust']['thrust'][1]))
    print('    @ %4.2f Isp\n' % opts['thrust']['Isp'])
    print('    @ %2.4e kg/s\n' % opts['thrust']['mdot'])
    print('\n. . . . . . . . . . . . . . . .\n\n')
    #%%
    return sorted_final_soln[0,:]