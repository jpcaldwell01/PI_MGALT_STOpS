# FORM: [J_fin] = MGALT_FBSM_costFun(J_init,pos,ang,vel_rad,vel_tan,tof,...
#       CONST,OPT,VAR)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function calculates the cost of a particular member by 
# |     looking at the results from the planetary transfers and comparing 
# |     them to the desired conditions
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -J_init             (1,1)       [float]         [unitless]
# |     	The cost of this member, denoted as 'f' in other functions
# |     -pos        	(2,transfers) 	[float]      	[DU]
# |         The radial positions of the spacecraft and planet (scplanet)
# |     -ang        	(2,transfers) 	[float]      	[deg]
# |         The angulat positions of the spacecraft and planet (scplanet)
# |     -vel_rad       	(2,transfers) 	[float]      	[DU/TU]
# |         The radial velocity of the spacecraft and planet (scplanet)
# |     -vel_tan       	(2,transfers) 	[float]      	[DU/TU]
# |         The angulat velocity of the spacecraft and planet (scplanet)
# |     -tof        	(2,transfers) 	[float]      	[TU]
# |         The times of flight for each transfer
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
# |     -J_fin              (1,1)       [float]         [unitless]
# |     	The cost of this member, denoted as 'f' in other functions
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------


import numpy as np

def MGALT_FBSM_costFun(J_init,pos,ang,vel_rad,vel_tan,tof,const,opts,var):

    ## Initials
    
    # Constants
    array_bodies = [3,4,5]
    J_trans = np.zeros((1,var['transfers']))
    J_trans_all = np.zeros((1,var['transfers']*5))
    
    
    
    ## Cost
    if var['transfers'] == 1:
        
        # ---------- Minimize Final Position Error ----------------------------
        JR = (pos[0] - pos[1])**2/(opts['cost']['tolR']**2)
        J_trans = J_trans + (opts['cost']['R']*JR)
    
        # ---------- Minimize Final Angular Displacement Error ----------------
        JT = (ang[0] - ang[1])**2/(opts['cost']['tolTheta']**2)
        J_trans = J_trans + (opts['cost']['Theta']*JT)
    
        # ---------- Minimize Final Radial Velocity Error ---------------------
        JU = (vel_rad[0] - vel_rad[1])**2/(opts['cost']['tolU']**2)
        J_trans = J_trans + (opts['cost']['U']*JU)
    
        # ---------- Minimize Final Tangential Velocity Error -----------------
        JV = (vel_tan[0] - vel_tan[1])**2/(opts['cost']['tolV']**2)
        J_trans = J_trans + (opts['cost']['V']*JV)
    
        # ---------- Restraint on Time of Flight ------------------------------
        TU = const['TU'][0]*86400                     # [sec/TU]
        tt_end = opts['thrust']['tt_end'][0]*86400*(1/TU)    # TU
        Jtt = (tof[0] + tof[1] - tt_end)**2/(opts['weighting']['W_tof_conv']**2)
        J_trans = J_trans + (opts['cost']['tt']*Jtt)
    
    else:
        for i4 in range(var['transfers']):
           
            # ---------- Minimize Final Position Error ----------------------------
            JR = (pos[i4,0] - pos[i4,1])**2/(opts['cost']['tolR']**2)
            J_trans[0,i4] = J_trans[0,i4] + (opts['cost']['R']*JR)
        
            # ---------- Minimize Final Angular Displacement Error ----------------
            JT = (ang[i4,0] - ang[i4,1])**2/(opts['cost']['tolTheta']**2)
            J_trans[0,i4] = J_trans[0,i4] + (opts['cost']['Theta']*JT)
        
            # ---------- Minimize Final Radial Velocity Error ---------------------
            JU = (vel_rad[i4,0] - vel_rad[i4,1])**2/(opts['cost']['tolU']**2)
            J_trans[0,i4] = J_trans[0,i4] + (opts['cost']['U']*JU)
        
            # ---------- Minimize Final Tangential Velocity Error -----------------
            JV = (vel_tan[i4,0] - vel_tan[i4,1])**2/(opts['cost']['tolV']**2)
            J_trans[0,i4] = J_trans[0,i4] + (opts['cost']['V']*JV)
        
            # ---------- Restraint on Time of Flight ------------------------------
            TU = const['TU'][0]*86400                     # [sec/TU]
            tt_end = opts['thrust']['tt_end'][i4]*86400*(1/TU)    # TU
            Jtt = (tof[i4,0] + tof[i4,1] - tt_end)**2/(opts['weighting']['W_tof_conv']**2)
            J_trans[0,i4] = J_trans[0,i4] + (opts['cost']['tt']*Jtt)
                    
    # Add all the costs together
    J_fin = np.sum(J_trans) + J_init
    
    return J_fin, J_trans
