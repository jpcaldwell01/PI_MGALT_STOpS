# FORM: [feasible,pert,transfer] = ...
#       MGALT_MBH_isFeasible(BOD,CONST,OPT,VAR,plot_vars,per_feas)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function looks at all of the transfers from the previous 
# |     iteration of the members within "MGALT_MBH_function". If the 
# |     transfers do not intersect the planet, then they are determined 
# |     to not be feasible. If the transfers do intersect the planet, 
# |     then the whole member is considered to be feasible.
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
# |     -plot_vars          (1,1)       [struct]     	[unitless]
# |         An object containing a lot of information about the 
# |         optimization parameters including: transfers(t and y ode 
# |         outputs), thrust values, thruster pointing angles, transfer 
# |         starting position, planet start/end locations for each 
# |         transfer, JD of each transfer, and tspans of each transfer
# |     -per_feas           (1,1)       [float]     	[unitless]
# |         The current feasibility percentage
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -feasible           (1,1)       [boolean]   	[unitless]
# |         The new perturbed member consisting of the original member and
# |         the perturbations added to it
# |     -pert  	(inputs.transfers,4)   	[float]       	[JD]
# |         The pert locations for this MBH struct
# |     -transfer        	(1,1)       [int]           [unitless]
# |         The current transfer number
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |     -For future improvement, this function should look at individual 
# |     segments of the transfer(s). If an individual segment is feasible, 
# |     then it should get flagged as a potential candidate and cycled 
# |     through into the cluster. The segments of the member which are 
# |     not feasible should then be altered until they are feasible.
# |
# |     -Some of the above infrastructure is setup with the pert and
# |     transfer returns, but they were not fully implimented.
# |
# |-----------------------------------------------------------------------

import numpy as np

def MGALT_MBH_isFeasible(bod,const,opts,opt_algo,var,plot_vars,per_feas):

    ## Initials
    high = np.zeros(np.shape(var['high']))
    low = np.zeros(np.shape(var['low']))

    for it in range(len(var['bin'])):
        high[it] = np.min([var['high'][it],var['highC'][it]])
        low[it] = np.max([var['low'][it],var['lowC'][it]])
        
    feasible = [False for x in range(var['transfers'])]
    transfer = [0 for x in range(var['transfers'])]
    count = 0
    
    # Get some constants
    AU = const['AU'][0]           #km/AU
    
    ## Disregard Departure and ToF Members
    
    # Compare against the times to see if random perturbations made the
    # transfer outside of the defined limit
    leave = np.zeros((var['transfers'],2))
    tof = np.zeros((var['transfers'],2))
    
    if opts['solver'] in ['LT_IN_FSM_2D','LT_DIR_FSM_2D']:
        '''
        # Run through all of the transfers
        for i1 = 1:VAR.transfers

            # Departure bounds
            leave(1,1) = plot_vars['JD'](1,1) < var['low'](1)
            leave(1,2) = plot_vars['JD'](1,1) > var['high'](1)

            # TOF Bounds
            tof(1,1) = plot_vars['JD'](1,2) < (var['low'](1)...
                + var['low'](end))
            tof(1,2) = plot_vars['JD'](1,2) > (var['high'](1)...
                + var['high'](end))

        end
        '''
    elif opts['solver'] in ['MGALT_IN_FBSM_2D','MGALT_DIR_FBSM_2D']:
        
        # Departure bounds
        leave[0,0] = plot_vars['JD'][0,0] < low[0]
        leave[0,1] = plot_vars['JD'][0,0] > high[0]
        
        if var['transfers'] == 1:
        
            tof[0,0] = plot_vars['JD'][0,1] < (low[0] + low[-1])
            tof[0,1] = plot_vars['JD'][0,1] > (high[0] + high[-1])
            
        else:
            
            tof[0,0] = plot_vars['JD'][0,1] < (low[0] + low[10])
            tof[0,1] = plot_vars['JD'][0,1] > (high[0] + high[10])
            
            # Run through all of the transfers
            for i1 in range(var['transfers']-2):
                leave[i1+1,0] = plot_vars['JD'][i1+1,0] < (low[i1*11+11])
                leave[i1+1,1] = plot_vars['JD'][i1+1,0] > (high[i1*11+11])
                tof[i1+1,0] = plot_vars['JD'][i1+1,1] < (low[i1*11+11] + low[i1*11+21])
                tof[i1+1,1] = plot_vars['JD'][i1+1,1] > (high[i1*11+11] + high[i1*11+21])
            
            leave[-1,0] = plot_vars['JD'][-1,0] < low[-8]
            leave[-1,1] = plot_vars['JD'][-1,0] > high[-8]
            tof[-1,0] = plot_vars['JD'][-1,1] < (low[-8] + low[-1])
            tof[-1,1] = plot_vars['JD'][-1,1] > (high[-8] + high[-1])
            
    else:
        
        raise Exception("invalid Solver Selection")
           
    pert = np.append(leave,tof)
    
    # The perturbations will break the ToF, exit
    if pert.any():
        return feasible,pert,transfer
    
    ## Compare Feasibilities
    
    if opts['solver'] in ['LT_IN_FSM_2D','LT_DIR_FSM_2D']:	# Feasibility is defined as intersecting with the planet's SOI
        '''
        # Check to see if intersects with SOI
        [does_intersect,pos_body,values] = ...
        SOIIntersectCheck(BOD,CONST,OPT,VAR,plot_vars,2)

        # If it does intersect with the SOI
        if does_intersect
            SOI_rad = CONST.(strcat(BOD.bodies{2},"_SOI"))/AU/2	# km/AU
            SOI_percent = ((SOI_rad*per_feas) + SOI_rad)         	# km/AU

            pos_sc = values(1,1:2)

            # Check to see if the terminal point is located within the
            # region defined as feasible
            if ( (norm(pos_sc) >= (norm(pos_body(1:2))-SOI_percent)) ...
                    && (norm(pos_sc) <= (norm(pos_body(1:2))+SOI_percent)) )

                feasible(1) = true
                transfer(1) = 1

            end

        end
        '''
    elif opts['solver'] == 'MGALT_IN_FBSM_2D':            	# Feasibility is defined as the patch points being within a certain tolerance of each other
        
        # Run through all of the transfers
        for i2 in range(var['transfers']):
            
            forward = plot_vars['transfers_fs'][:,int((i2)*(np.shape(plot_vars['transfers_fs'])[1]/var['transfers'])):int((i2+1)*(np.shape(plot_vars['transfers_fs'])[1]/var['transfers']))]
            backward = plot_vars['transfers_bs'][:,int((i2)*(np.shape(plot_vars['transfers_fs'])[1]/var['transfers'])):int((i2+1)*(np.shape(plot_vars['transfers_fs'])[1]/var['transfers']))]
            
            # I need to see if the match points are feasible.
            # Dr. Englander recommended using 1e-5 AU for the tolerance
            # point
            tol = opt_algo['feas_tol']
            check = np.zeros((3,))
            
            # Get the position norms, vel rad, and vel tan
            # norm_pos_fs = np.linalg.norm([forward[-1,4],forward[-1,5]])
            # norm_pos_bs = np.linalg.norm([backward[-1,4],backward[-1,5]])
            norm_pos_fs = 1*forward[-1,4]
            norm_pos_bs = 1*backward[-1,4]
            
            vel_rad_fs = forward[-1,6]
            vel_rad_bs = backward[-1,6]
            
            vel_tan_fs = forward[-1,7]
            vel_tan_bs = backward[-1,7]
            
            # Check feasibility of position X and Y
            if (abs(norm_pos_fs-norm_pos_bs) <= (tol*per_feas)):
                
                # Make sure in the same quadrant
                if (np.sign(forward[-1,4]) == np.sign(backward[-1,4])) and (np.sign(forward[-1,5]) == np.sign(backward[-1,5])):
                    check[0] = 1
            
            # Check feasibility of radial velocity
            if (abs(vel_rad_fs-vel_rad_bs) <= (tol*per_feas)):
                check[1] = 1
            
            # Check feasibility of tangential velocity
            if (abs(vel_tan_fs-vel_tan_bs) <= (tol*per_feas)):
                check[2] = 1
            
            # If all of the checks pass
            if check.all():
                feasible[count] = True
                transfer[count] = i2
            else:
                feasible[count] = False
                transfer[count] = i2
            
            count = count+1
            
    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':            	# Feasibility is defined as the patch points being within a certain tolerance of each other
        
        # Run through all of the transfers
        for i2 in range(var['transfers']):
            
            forward = plot_vars['transfers_fs'][:,int((i2)*(np.shape(plot_vars['transfers_fs'])[1]/var['transfers'])):int((i2+1)*(np.shape(plot_vars['transfers_fs'])[1]/var['transfers']))]
            backward = plot_vars['transfers_bs'][:,int((i2)*(np.shape(plot_vars['transfers_fs'])[1]/var['transfers'])):int((i2+1)*(np.shape(plot_vars['transfers_fs'])[1]/var['transfers']))]
            
            # I need to see if the match points are feasible.
            # Dr. Englander recommended using 1e-5 AU for the tolerance
            # point
            tol = opt_algo['feas_tol']
            check = np.zeros((3,))
            
            # Get the position norms, vel rad, and vel tan
            norm_pos_fs = 1*forward[-1,1]
            norm_pos_bs = 1*backward[-1,1]
            
            vel_rad_fs = 1*forward[-1,3]
            vel_rad_bs = 1*backward[-1,3]
            
            vel_tan_fs = 1*forward[-1,4]
            vel_tan_bs = 1*backward[-1,4]
            
            # Check feasibility of position X and Y
            if abs(norm_pos_fs-norm_pos_bs) <= (tol*per_feas):
                
                # Make sure in the same quadrant
                if np.sign(forward[-1,1]) == np.sign(backward[-1,1]) and np.sign(forward[-1,2]) == np.sign(backward[-1,2]):
                    check[0] = 1
            
            # Check feasibility of radial velocity
            if abs(vel_rad_fs-vel_rad_bs) <= (tol*per_feas):
                check[0] = 1
            
            # Check feasibility of tangential velocity
            if abs(vel_tan_fs-vel_tan_bs) <= (tol*per_feas):
                check[2] = 1
            
            # If all of the checks pass
            if check.all():
                feasible[count] = True
                transfer[count] = i2
            else:
                feasible[count] = False
                transfer[count] = i2
                
            count = count+1
            
    else:
        raise Exception('Invalid Solver Selection')
    
    return feasible,pert,transfer