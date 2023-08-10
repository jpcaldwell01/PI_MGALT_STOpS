# FORM: [R_init,V_init,pos_rad,pos_ang,vel_rad,vel_tan] = ...
#       MGALT_conditionsInit(JD,BOD,CONST,OPT,VAR,planet_parse)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function gets a lot of different parameters for the 
# |     departure planet. These include the R and V of the departure at 
# |     the first Julian Date, as well as the radial position and 
# |     velocity and the angular location and velocity
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -JD                 (1,1)       [float]         [JD]
# |         Julian day for calculation at departure
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
# |     -planet_parse       (3,1)       [int]        [unitless]
# |         An array containing the index numbers for the bodies
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -R_init             (3,1)       [float]         [AU]
# |         Radius vector
# |     -V_init             (3,1)       [float]         [AU/TU]
# |         Velocity vector at JD_1 and JD_2
# |     -pos_rad        	(1,1)   	[float]         [DU]
# |         The planet radial position
# |     -pos_ang        	(1,1)   	[float]         [deg]
# |         The planet angular position
# |     -vel_rad        	(1,1)   	[float]         [DU/TU]
# |         The planet radial velocity
# |     -vel_tan        	(1,1)   	[float]         [DU/TU]
# |         The planet tangential velocity
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

from Other_Functions.MGALT_stateBodies import MGALT_stateBodies
import math
import numpy as np

def MGALT_conditionsInit(JD,bod,const,opts,var,planet_parse):

    ## Solve initial spacecraft conditions
    
    # Initials
    AU = const['AU'][0]        	# [km/AU]
    TU = const['TU'][0]*86400	# [sec/TU]
    
    
    # Get R and V for the specified planet at the specified JD
    [R_init,V_init] = MGALT_stateBodies(JD,bod,const,opts,var,planet_parse)
    
    # For only the first planet, get the initial R and V vectors
    R0 = R_init*(1/AU)             # DU
    V0 = V_init*(TU/AU)            # TU
    
    # Convert into polar coords
    pos_rad = np.linalg.norm(R0)             # DU
    pos_ang = math.degrees(math.atan2(R0[1],R0[0]))	# (deg) ccw from positive X
    
    # Account for tangent region
    if pos_ang < 0:
        pos_ang = 360 + pos_ang
    
    
    # Get the coords in LVLH
    R0_lvlh = R0/pos_rad
    W0_lvlh = np.cross(R0,V0)/np.linalg.norm(np.cross(R0,V0))
    S0_lvlh = np.cross(W0_lvlh,R0_lvlh)
    T0_ECI_lvlh = np.stack((R0_lvlh, S0_lvlh, W0_lvlh),axis=1)
    V0_lvlh = np.dot(T0_ECI_lvlh.transpose(),V0)     # DU/TU
    
    # Radial and Tangential Velocity for the solver
    vel_rad = V0_lvlh[0]           # DU/TU, parallel to R vector
    vel_tan = V0_lvlh[1]           # DU/TU, perpendicular to R vector
    
    return R_init,V_init,pos_rad,pos_ang,vel_rad,vel_tan
