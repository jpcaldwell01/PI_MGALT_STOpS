# FORM: [J,plot_vars] = MGALT_IN_FBSM_2D(member,BOD,CONST,OPT,VAR)
#
# |-----------------------------------------------------------------------
# | NOTES:
# |     -Cost function solver for use with MGALT STOpS
# |
# |     -This is the INDIRECT method for representing a low-thrust 
# |     orbital tragectory, taken from Conway
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -member             (1,Nvar)    [float]         [unitless]
# |         A single member of a population input into the function
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
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -J                  (1,1)       [float]         [unitless]
# |     	The cost of this member, denoted as 'f' in other functions
# |     -plot_vars          (1,1)       [struct]     	[unitless]
# |         An object containing a lot of information about the 
# |         optimization parameters including: transfers(t and y ode 
# |         outputs), thrust values, thruster pointing angles, transfer 
# |         starting position, planet start/end locations for each 
# |         transfer, JD of each transfer, and tspans of each transfer
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |     -References
# |         B. A. Conway, Spacecraft Trajectory Optimization. Cambridge University Press, 2010.
# |
# |-----------------------------------------------------------------------

# from getConstants import getConstants
import numpy as np
from Solvers.Transfer_Conditions.MGALT_conditionsInit import MGALT_conditionsInit
from Solvers.Indirect.MGALT_IN_FBSM_2D_EOM import MGALT_IN_FBSM_2D_EOM
from Solvers.Transfer_Conditions.MGALT_conditionsTransFBSM import MGALT_conditionsTransFBSM
from scipy.integrate import solve_ivp
import math
from Solvers.Transfer_Conditions.MGALT_conditionsGravityAssist import MGALT_conditionsGravityAssist
from Solvers.Transfer_Conditions.MGALT_convertHelioLVLH import MGALT_convertHelioLVLH
from Solvers.Cost_Functions.MGALT_FBSM_costFun import MGALT_FBSM_costFun

def MGALT_IN_FBSM_2D(member,bod,const,opts,var):

    ## Setup
    
    # Canonical Units
    AU = const['AU'][0]                  # [km/AU]
    TU = const['TU'][0]*86400            # [sec/TU]
    mew_sun = const['Sun_mu'][0]         # [km^3/s^2]
    mew_sc = mew_sun*(TU**2/AU**3)  	# [DU^3/TU^2] Mew for the spacecraft
    
    # Total Cost
    J = 0
    
    # ODE Parameters
    # options = odeset('AbsTol',1e-6,'RelTol',1e-6)
    tspan_divider = 400
    
    # Other
    transfers = var['transfers']
    
    ## Pre-Allocation
    
    plot_vars = {}
    plot_vars['Y0_fs'] = np.zeros((8,transfers))
    plot_vars['Y0_bs'] = np.zeros((8,transfers))
    plot_vars['transfers_fs'] = np.zeros((tspan_divider,9*transfers))
    plot_vars['transfers_bs'] = np.zeros((tspan_divider,9*transfers))
    plot_vars['transfers'] = np.zeros((2*tspan_divider,9*transfers))
    tspan_fs = np.zeros((1,tspan_divider))
    tspan_bs = np.zeros((1,tspan_divider))
    
    # Orbit variables
    tspan = np.zeros((transfers,2))
    JD = np.zeros((transfers,2))
    
    # Cost function variables
    pos_rad_sc_fs = np.zeros((1,transfers))
    pos_rad_sc_bs = np.zeros((1,transfers))
    pos_ang_sc_fs = np.zeros((1,transfers))
    pos_ang_sc_bs = np.zeros((1,transfers))
    vel_rad_sc_fs = np.zeros((1,transfers))
    vel_rad_sc_bs = np.zeros((1,transfers))
    vel_tan_sc_fs = np.zeros((1,transfers))
    vel_tan_sc_bs = np.zeros((1,transfers))
    transfer_time_sc_fs = np.zeros((1,transfers))
    transfer_time_sc_bs = np.zeros((1,transfers))
    
    
    # Transfer variables
    planet_departure = np.zeros((6,1))
    planet_trans = np.zeros((6,transfers-1))
    planet_target = np.zeros((6,1))
    
    
    
    ## Perform Transfer(s)
    
    if transfers == 1:        # Going from planet A to B

        # JD for the bodies during transfer segments
        JD = np.array([member[0], member[0]+member[-1]])     # Start and end Julian Day
        JD.shape = (1,2)

        # tspan for all transfer segments
        tspan = np.array([value*(1/TU) for value in [0, member[-1]*86400]])        # TU
        tspan.shape = (1,2)

        # ********** BODY CONDITIONS **********
        # Get departure body initial conditions
        [planet_R_dep,planet_V_dep,sc_dep_pos_rad,sc_dep_pos_ang,sc_dep_vel_rad,sc_dep_vel_tan] = MGALT_conditionsInit(JD[0,0],bod,const,opts,var,[0,1,2])
        
        # Get target body final conditions
        [planet_R_tar,planet_V_tar,sc_tar_pos_rad,sc_tar_pos_ang,sc_tar_vel_rad,sc_tar_vel_tan] = MGALT_conditionsInit(JD[-1,-1],bod,const,opts,var,[3,4,5])
        
        
        # ********** FORWARD SHOOTING **********
        # Check to see if any additional dV from launch vehicle
        if 'launch_dV_rad' in opts['thrust']:
            sc_dep_vel_rad = sc_dep_vel_rad + (opts['thrust']['launch_dV_rad']*(TU/AU))
        
        if 'launch_dV_tan' in opts['thrust']:
            sc_dep_vel_tan = sc_dep_vel_tan + (opts['thrust']['launch_dV_tan']*(TU/AU))
        
        
        # Define the coefficients for forward shooting
        lam1_fs = member[1]
        lam2_fs = member[2]
        lam3_fs = member[3]
        mass_fs_Y0 = opts['thrust']['m0']     # inital mass (kg)

        # Define the state varaible for forward shooting
        Y0_fs = [lam1_fs,lam2_fs,lam3_fs,sc_dep_pos_rad,sc_dep_pos_ang,sc_dep_vel_rad,sc_dep_vel_tan,mass_fs_Y0]
        tspan_fs = np.linspace(tspan[0,0],tspan[0,-1]/2,tspan_divider)


        # ODE45 to solve for forward segment
        # [ttot_fs,Ytot_fs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
        #     tspan_fs,Y0_fs,OPT.ode,CONST,OPT,mew_sc)
        
        try:
            solution_fs = solve_ivp(MGALT_IN_FBSM_2D_EOM,[tspan_fs[0],tspan_fs[-1]],Y0_fs,t_eval=tspan_fs,args=(const,opts,mew_sc),rtol=1e-8)
            Ytot_fs = solution_fs.y.T
            ttot_fs = solution_fs.t.T
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            
            return J,plot_vars,J_trans
        
        # End Conditions
        pos_rad_sc_fs      	= Ytot_fs[-1,3]
        pos_ang_sc_fs     	= (Ytot_fs[-1,4]/360 - math.floor(Ytot_fs[-1,4]/360))*360   # Accounts for being larger than 360
        vel_rad_sc_fs    	    = Ytot_fs[-1,5]
        vel_tan_sc_fs     	= Ytot_fs[-1,6]
        transfer_time_sc_fs = tspan_fs[0]
       
        
        # ********** BACKWARDS SHOOTING **********
        # Define the coefficients for backward shooting
        lam1_bs = member[-4]
        lam2_bs = member[-5]
        lam3_bs = member[-6]
        
        # Assume a constant mass loss for the duration of the transfer
        mass_bs_Y0 = mass_fs_Y0 - 2*(Ytot_fs[0,-1]-Ytot_fs[-1,-1])

        # Define the state varaible for backwards shooting
        Y0_bs = [lam1_bs,lam2_bs,lam3_bs,sc_tar_pos_rad,sc_tar_pos_ang,sc_tar_vel_rad,sc_tar_vel_tan,mass_bs_Y0]
        tspan_bs = np.linspace(tspan[0,-1]/2,tspan[0,0],tspan_divider)

        # ODE45 to solve for backwards segment
        # [ttot_bs,Ytot_bs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
            # tspan_bs,Y0_bs,OPT.ode,CONST,OPT,mew_sc)
        
        try:
            solution_bs = solve_ivp(MGALT_IN_FBSM_2D_EOM,[tspan_bs[0],tspan_bs[-1]],Y0_bs,t_eval=tspan_bs,args=(const,opts,mew_sc),rtol=1e-8)
            Ytot_bs = solution_bs.y.T
            ttot_bs = solution_bs.t.T
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            return J,plot_vars,J_trans
        
        # End Conditions
        pos_rad_sc_bs      	= Ytot_bs[-1,3]
        pos_ang_sc_bs     	= (Ytot_bs[-1,4]/360 - math.floor(Ytot_bs[-1,4]/360))*360   # Accounts for being larger than 360
        vel_rad_sc_bs    	    = Ytot_bs[-1,5]
        vel_tan_sc_bs     	= Ytot_bs[-1,6]
        transfer_time_sc_bs = tspan_bs[0]

        
        # ********** Append to spacecraft plotting variables **********
        # planet_departure = [planet_R_depplanet_V_dep]
        # planet_target = [planet_R_tarplanet_V_tar]
        # plot_vars.Y0_fs{1,1} = Y0_fs
        # plot_vars.Y0_bs{1,1} = Y0_bs
        # plot_vars.transfers_fs{1,1} = [ttot_fs,Ytot_fs]
        # plot_vars.transfers_bs{1,1} = [ttot_bs,Ytot_bs]
        # plot_vars.transfers{1,1} = [ttot_fs,Ytot_fs...
        #                             flipud(ttot_bs)+ttot_fs(end),flipud(Ytot_bs)]
        try:
            planet_departure = np.append(planet_R_dep,planet_V_dep)
            planet_target = np.append(planet_R_tar,planet_V_tar)    # Plot vars
            plot_vars['Y0_fs'] = Y0_fs
            plot_vars['Y0_bs'] = Y0_bs
            plot_vars['transfers_fs'] = np.hstack((tspan_fs.T[:,None],Ytot_fs))
            plot_vars['transfers_bs'] = np.hstack((tspan_bs.T[:,None],Ytot_bs))
            plot_vars['transfers'] = np.vstack( (np.hstack((tspan_fs.T[:,None],Ytot_fs)),np.hstack((np.flipud(tspan_bs[:,None])+tspan_fs[-1],np.flipud(Ytot_bs)))))
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            
            return J,plot_vars,J_trans
        
    else:       # Going from planet A to N via B, C, D, ...
        
        # Misc
        cost_count = 0

        # For indexing the positions
        array_bodies = [3,4,5]
        array_member = [8,9]
        
        # JD for the bodies during transfer segments
        start = member[0]
        tof = member[10]
        JD[0] = [start, start+tof]
        for i1 in range(1,transfers-1):
            start = member[i1*11]
            tof = member[i1*11+10]
            JD[i1] = [start, start+tof]                  # Start and end Julian Day
        
        start = member[-8]
        tof = member[-1]
        JD[-1] = [start, start+tof]
    
        
        # tspan for all transfer segments
        tspan[0] = [value*(1/TU) for value in [0, member[10]*86400]]
        for i2 in range(1,transfers-1):
            tspan[i2] = [value*(1/TU) for value in [0, member[i2*11+10]*86400]]#TU
        
        tspan[-1] = [value*(1/TU) for value in [0, member[-1]*86400]]
        
        ## ********** DEPARTURE -> TRANSFER 1 **********
        
        # Departure body initial conditions
        [planet_R_dep,planet_V_dep,sc_dep_pos_rad,sc_dep_pos_ang,sc_dep_vel_rad,sc_dep_vel_tan] = MGALT_conditionsInit(JD[0,0],bod,const,opts,var,[0,1,2])
        
        # Check to see if any additional dV from launch vehicle
        if 'launch_dV_rad' in opts['thrust'].keys():
            sc_dep_vel_rad = sc_dep_vel_rad + (opts['thrust']['launch_dV_rad']*(TU/AU))
        if 'launch_dV_tan' in opts['thrust'].keys():
            sc_dep_vel_tan = sc_dep_vel_tan + (opts['thrust']['launch_dV_tan']*(TU/AU))
        
        # Define the state variable for forward shooting
        lam1_fs = member[1]
        lam2_fs = member[2]
        lam3_fs = member[3]
        mass_fs_Y0 = opts['thrust']['m0']     # inital mass (kg)
        
        Y0_fs = [lam1_fs,lam2_fs,lam3_fs,sc_dep_pos_rad,sc_dep_pos_ang,sc_dep_vel_rad,sc_dep_vel_tan,mass_fs_Y0]
        tspan_fs = np.linspace(tspan[0,0],tspan[0,-1]/2,tspan_divider)
        
        # ODE45 to solve for forward segment
        # [ttot_fs,Ytot_fs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
            # tspan_fs,Y0_fs,OPT.ode,CONST,OPT,mew_sc)
        try:
            solution_fs = solve_ivp(MGALT_IN_FBSM_2D_EOM,[tspan_fs[0],tspan_fs[-1]],Y0_fs,t_eval=tspan_fs,args=(const,opts,mew_sc),rtol=1e-8)
            Ytot_fs = solution_fs.y.T
            ttot_fs = solution_fs.t.T
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            
            return J,plot_vars,J_trans
        
        # FS end Conditions
        pos_rad_sc_fs[0,cost_count]     	= Ytot_fs[-1,3]
        pos_ang_sc_fs[0,cost_count]    	= (Ytot_fs[-1,4]/360 - math.floor(Ytot_fs[-1,4]/360))*360   # Accounts for being larger than 360
        vel_rad_sc_fs[0,cost_count]   	= Ytot_fs[-1,5]
        vel_tan_sc_fs[0,cost_count]     	= Ytot_fs[-1,6]
        mass_sc_fs_end                  = Ytot_fs[-1,7]
        transfer_time_sc_fs[0,cost_count]	= ttot_fs[-1]
        
        # Get transfer body position and s/c positions
        # [planet_R_trans,planet_V_trans,sc_trans_pos_rad,sc_trans_pos_ang,sc_trans_vel_rad,sc_trans_vel_tan,control,sc_vel_helio_enter]
        [planet_R_trans,planet_V_trans,sc_trans_pos_rad,sc_trans_pos_ang,sc_trans_vel_rad,sc_trans_vel_tan,control,sc_vel_helio_enter] = MGALT_conditionsTransFBSM(JD[0,1],bod,const,opts,var,tspan[0],mass_sc_fs_end,[1,1],member[array_member].tolist() + [opts['weighting']['control_v']],sc_dep_pos_rad,array_bodies)

        # Change the control numbers
        member[array_member] = control[0:2]
        
        # Define the state varaible for backwards shooting
        lam1_bs = member[4]
        lam2_bs = member[5]
        lam3_bs = member[6]
        mass_bs_Y0 = mass_fs_Y0 - 2*(Ytot_fs[0,-1]-Ytot_fs[-1,-1])
        
        Y0_bs = [lam1_bs,lam2_bs,lam3_bs,sc_trans_pos_rad,sc_trans_pos_ang,sc_trans_vel_rad,sc_trans_vel_tan,mass_bs_Y0]
        tspan_bs = np.linspace(tspan[0,-1]/2,tspan[0,0],tspan_divider)
        
        # ODE45 to solve for backwards segment
        try:
            solution_bs = solve_ivp(MGALT_IN_FBSM_2D_EOM,[tspan_bs[0],tspan_bs[-1]],Y0_bs,t_eval=tspan_bs,args=(const,opts,mew_sc),rtol=1e-8)
            Ytot_bs = solution_bs.y.T
            ttot_bs = solution_bs.t.T
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            return J,plot_vars,J_trans
        
        # BS end Conditions
        pos_rad_sc_bs[0,cost_count]      	= Ytot_bs[-1,3]
        pos_ang_sc_bs[0,cost_count]     	= (Ytot_bs[-1,4]/360 - math.floor(Ytot_bs[-1,4]/360))*360   # Accounts for being larger than 360
        vel_rad_sc_bs[0,cost_count]    	= Ytot_bs[-1,5]
        vel_tan_sc_bs[0,cost_count]     	= Ytot_bs[-1,6]
        transfer_time_sc_bs[0,cost_count] = tspan_bs[0]
        

        # Count up the array bodies
        array_bodies = [val+3 for val in array_bodies]
        array_member = [val+11 for val in array_member]
        cost_count = cost_count+1
        
        # Append to spacecraft plotting variables
        try:
            planet_departure = np.append(planet_R_dep,planet_V_dep)
            planet_trans[:,0] = np.append(planet_R_trans,planet_V_trans)    # Plot vars
            plot_vars['Y0_fs'][:,0] = Y0_fs
            plot_vars['Y0_bs'][:,0] = Y0_bs
            plot_vars['transfers_fs'][:,0:9] = np.hstack((tspan_fs.T[:,None],Ytot_fs))
            plot_vars['transfers_bs'][:,0:9] = np.hstack((tspan_bs.T[:,None],Ytot_bs))
            plot_vars['transfers'][:,0:9] = np.vstack( (np.hstack((tspan_fs.T[:,None],Ytot_fs)),np.hstack((np.flipud(tspan_bs[:,None])+tspan_fs[-1],np.flipud(Ytot_bs)))))
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            
            return J,plot_vars,J_trans
        

        ##  ********** TRANSFER 1 -> TRANSFER n **********
        
        for i3 in range(1,transfers-1):
            
            # Calculate the FS departure position and vel
            rp_coe = member[array_member[0]-12]        # Number between 0 and 1 for flyby altitude
            
            # Get the gravity assist condition
            
            sc_vel_helio_exit = MGALT_conditionsGravityAssist(bod,const,planet_R_trans,planet_V_trans,sc_vel_helio_enter,rp_coe,i3)[0]
            
            # Convert flyby LVLH
            sc_dep_vel = MGALT_convertHelioLVLH(const,planet_R_trans,sc_vel_helio_exit)[1]

            # Define the state varaible for forward shooting
            lam1_fs = member[-7]
            lam2_fs = member[-6]
            lam3_fs = member[-5]
            mass_fs_Y0 = mass_bs_Y0    # kg

            Y0_fs = [lam1_fs,lam2_fs,lam3_fs,sc_trans_pos_rad,sc_trans_pos_ang,sc_dep_vel[0],sc_dep_vel[1],mass_fs_Y0]
            tspan_fs = np.linspace(tspan[i3,0],tspan[i3,-1]/2,tspan_divider)

            # ODE45 to solve for forward segment
            try:
                solution_fs = solve_ivp(MGALT_IN_FBSM_2D_EOM,[tspan_fs[0],tspan_fs[-1]],Y0_fs,t_eval=tspan_fs,args=(const,opts,mew_sc),rtol=1e-8)
                Ytot_fs = solution_fs.y.T
                ttot_fs = solution_fs.t.T
            except:
                J = 9999999999999999
                J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
                
                return J,plot_vars,J_trans
            
            # FS end Conditions
            pos_rad_sc_fs[0,cost_count]     	= Ytot_fs[-1,3]
            pos_ang_sc_fs[0,cost_count]    	= (Ytot_fs[-1,4]/360 - math.floor(Ytot_fs[-1,4]/360))*360   # Accounts for being larger than 360
            vel_rad_sc_fs[0,cost_count]   	= Ytot_fs[-1,5]
            vel_tan_sc_fs[0,cost_count]     	= Ytot_fs[-1,6]
            transfer_time_sc_fs[0,cost_count]	= tspan_fs[-1]
            
            
            # Get transfer body position and s/c positions                           
            [planet_R_trans,planet_V_trans,sc_trans_pos_rad,sc_trans_pos_ang,sc_trans_vel_rad,sc_trans_vel_tan,control,sc_vel_helio_enter] = MGALT_conditionsTransFBSM(JD[i3,1],bod,const,opts,var,tspan[i3],mass_sc_fs_end,[1,1],member[array_member].tolist() + [opts['weighting']['control_v']],sc_dep_pos_rad,array_bodies)

            
            # Change the control numbers
            member[array_member] = [control[0],control[1]]

            # Define the state varaible for backwards shooting
            lam1_bs = member[4]
            lam2_bs = member[5]
            lam3_bs = member[6]
            mass_bs_Y0 = mass_fs_Y0 - 2*(Ytot_fs[0,-1]-Ytot_fs[-1,-1])

            Y0_bs = [lam1_bs,lam2_bs,lam3_bs,sc_trans_pos_rad,sc_trans_pos_ang,sc_trans_vel_rad,sc_trans_vel_tan,mass_bs_Y0]
            tspan_bs = np.linspace(tspan[i3,-1]/2,tspan[i3,0],tspan_divider)
            
            # ODE45 to solve for backwards segment
            try:
                solution_bs = solve_ivp(MGALT_IN_FBSM_2D_EOM,[tspan_bs[0],tspan_bs[-1]],Y0_bs,t_eval=tspan_bs,args=(const,opts,mew_sc),rtol=1e-8)
                Ytot_bs = solution_bs.y.T
                ttot_bs = solution_bs.t.T
            except:
                J = 9999999999999999
                J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
                
                return J,plot_vars,J_trans
            
            # BS end Conditions
            pos_rad_sc_bs[0,cost_count]      	= Ytot_bs[-1,3]
            pos_ang_sc_bs[0,cost_count]     	= (Ytot_bs[-1,4]/360 - math.floor(Ytot_bs[-1,4]/360))*360   # Accounts for being larger than 360
            vel_rad_sc_bs[0,cost_count]    	= Ytot_bs[-1,5]
            vel_tan_sc_bs[0,cost_count]     	= Ytot_bs[-1,6]
            transfer_time_sc_bs[0,cost_count] = tspan_bs[0]

            # Count up the array bodies
            array_bodies = [val+3 for val in array_bodies]
            array_member = [val+11 for val in array_member]
            cost_count = cost_count+1
            
            # Append to spacecraft plotting variables
            try:
                planet_trans[:,i3] = np.append(planet_R_trans,planet_V_trans)    # Plot vars
                plot_vars['Y0_fs'][:,i3] = Y0_fs
                plot_vars['Y0_bs'][:,i3] = Y0_bs
                plot_vars['transfers_fs'][:,i3*9:(i3+1)*9] = np.hstack((tspan_fs.T[:,None],Ytot_fs))
                plot_vars['transfers_bs'][:,i3*9:(i3+1)*9] = np.hstack((tspan_bs.T[:,None],Ytot_bs))
                plot_vars['transfers'][:,i3*9:(i3+1)*9] = np.vstack( (np.hstack((tspan_fs.T[:,None],Ytot_fs)),np.hstack((np.flipud(tspan_bs[:,None])+tspan_fs[-1],np.flipud(Ytot_bs)))))
            except:
                J = 9999999999999999
                J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
                
                return J,plot_vars,J_trans
            
        ## ********** TRANSFER n -> TARGET **********
        
        # Calculate the FS departure position and vel
        rp_coe = member[array_member[0]-12]        # Number between 0 and 1 for flyby altitude
        
        # Get the gravity assist conditions
        sc_vel_helio_exit = MGALT_conditionsGravityAssist(bod,const,planet_R_trans,planet_V_trans,sc_vel_helio_enter,rp_coe,transfers-1)[0]
        
        sc_dep_vel = MGALT_convertHelioLVLH(const,planet_R_trans,sc_vel_helio_exit)[1]
        
        # Define the state variable for forward shooting
        lam1_fs = member[-7]
        lam2_fs = member[-6]
        lam3_fs = member[-5]
        mass_fs_Y0 = mass_bs_Y0    # kg
        
        Y0_fs = [lam1_fs,lam2_fs,lam3_fs,sc_trans_pos_rad,sc_trans_pos_ang,sc_dep_vel[0],sc_dep_vel[1],mass_fs_Y0]
        tspan_fs = np.linspace(tspan[-1,0],tspan[-1,-1]/2,tspan_divider)

        # ODE45 to solve for forward segment
        try:
            solution_fs = solve_ivp(MGALT_IN_FBSM_2D_EOM,[tspan_fs[0],tspan_fs[-1]],Y0_fs,t_eval=tspan_fs,args=(const,opts,mew_sc),rtol=1e-8)
            Ytot_fs = solution_fs.y.T
            ttot_fs = solution_fs.t.T
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            
            return J,plot_vars,J_trans
        
        # FS end Conditions
        pos_rad_sc_fs[0,cost_count]     	= Ytot_fs[-1,3]
        pos_ang_sc_fs[0,cost_count]    	= (Ytot_fs[-1,4]/360 - math.floor(Ytot_fs[-1,4]/360))*360   # Accounts for being larger than 360
        vel_rad_sc_fs[0,cost_count]   	= Ytot_fs[-1,5]
        vel_tan_sc_fs[0,cost_count]     	= Ytot_fs[-1,6]
        # mass_sc_fs_end                  = Ytot_fs[-1,7]
        transfer_time_sc_fs[0,cost_count]	= tspan_fs[-1]
        
        # Get target body final conditions
        [planet_R_tar,planet_V_tar,sc_tar_pos_rad,sc_tar_pos_ang,sc_tar_vel_rad,sc_tar_vel_tan] = MGALT_conditionsInit(JD[-1,-1],bod,const,opts,var,array_bodies)
        
        
        # Define the coefficients for backward shooting
        lam1_bs = member[-4]
        lam2_bs = member[-3]
        lam3_bs = member[-2]
        mass_bs_Y0 = mass_fs_Y0 - 2*(Ytot_fs[0,-1]-Ytot_fs[-1,-1])

        # Define the state varaible for backwards shooting
        Y0_bs = [lam1_bs,lam2_bs,lam3_bs,sc_tar_pos_rad,sc_tar_pos_ang,sc_tar_vel_rad,sc_tar_vel_tan,mass_bs_Y0]
        tspan_bs = np.linspace(tspan[-1,-1]/2,tspan[-1,0],tspan_divider)

        # ODE45 to solve for backwards segment
        try:
            solution_bs = solve_ivp(MGALT_IN_FBSM_2D_EOM,[tspan_bs[0],tspan_bs[-1]],Y0_bs,t_eval=tspan_bs,args=(const,opts,mew_sc),rtol=1e-8)
            Ytot_bs = solution_bs.y.T
            ttot_bs = solution_bs.t.T
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            
            return J,plot_vars,J_trans
        
        # End Conditions
        
        pos_rad_sc_bs[0,-1]      	= Ytot_bs[-1,3]
        pos_ang_sc_bs[0,-1]     	= (Ytot_bs[-1,4]/360 - math.floor(Ytot_bs[-1,4]/360))*360   # Accounts for being larger than 360
        vel_rad_sc_bs[0,-1]    	    = Ytot_bs[-1,5]
        vel_tan_sc_bs[0,-1]     	= Ytot_bs[-1,6]
        transfer_time_sc_bs[0,-1] = tspan_bs[0]
        
        # Plot Vars
        try:
            planet_target = np.append(planet_R_tar,planet_V_tar)
            plot_vars['Y0_fs'][:,-1] = Y0_fs
            plot_vars['Y0_bs'][:,-1] = Y0_bs
            plot_vars['transfers_fs'][:,-9:] = np.hstack((tspan_fs.T[:,None],Ytot_fs))
            plot_vars['transfers_bs'][:,-9:] = np.hstack((tspan_bs.T[:,None],Ytot_bs))
            plot_vars['transfers'][:,-9:] = np.vstack( (np.hstack((tspan_fs.T[:,None],Ytot_fs)),np.hstack((np.flipud(tspan_bs[:,None])+tspan_fs[-1],np.flipud(Ytot_bs)))))
        except:
            J = 9999999999999999
            J_trans = np.array([J/(len(bod['bodies'])-1) for i in range(len(bod['bodies'])-1)],ndmin=2)            
            
            return J,plot_vars,J_trans
    
    
    ## Plotting Vars and Difference
    
    # Plotting variables
    plot_vars['planetary_conditions'] = np.hstack((planet_departure[:,None], planet_trans, planet_target[:,None]))
    plot_vars['JD'] = JD          # variables necessary to plot the thrust vectors
    plot_vars['tspan'] = tspan  	# time spans for the orbits
    
    # plotOrbits(BOD,CONST,OPT,VAR,plot_vars)    
    
    
    ## Calculating Cost
    
    # Get the total cost function for the population member
    J,J_trans = MGALT_FBSM_costFun(J, np.hstack((pos_rad_sc_fs.T, pos_rad_sc_bs.T)),np.hstack((pos_ang_sc_fs.T, pos_ang_sc_bs.T)),np.hstack((vel_rad_sc_fs.T, vel_rad_sc_bs.T)),np.hstack((vel_tan_sc_fs.T, vel_tan_sc_bs.T)),np.hstack((transfer_time_sc_fs.T, transfer_time_sc_bs.T)),const,opts,var)
    
    return J,plot_vars,J_trans