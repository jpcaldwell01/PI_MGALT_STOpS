# FORM: [dY] = MGALT_IN_FBSM_2D_EOM(t,Y,CONST,OPT,mew_sc)
#
# |-----------------------------------------------------------------------
# | NOTES:
# |     -Equations of motion for Conway variable construction for 
# |     low-thrust trajectory optimization. 2D assumption
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -t                  (1,1)       [float]         [unitless]
# |         Time
# |     -Y                  (8,1)       [float]         [unitless]
# |         State vector, see MISC for unpacking
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
# |     -mew_sc           	(1,1)     	[float]         [DU^3/TU^2]
# |         Gravitaional parameter of the spacecraft
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -dY                 (8,1)       [float]         [unitless]
# |         Derivative of state vector
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |     -Unpacking
# |         Y(0) = lamda1      	(unitless) costate variable 1
# |         Y(1) = lamda2      	(unitless) costate variable 2
# |         Y(2) = lamda3   	(unitless) costate variable 3
# |         Y(3) = pos_rad      (DU) radial position
# |         Y(4) = pos_ang      (deg) angular position
# |         Y(5) = vel_rad      (DU/TU) radial velocity
# |         Y(6) = vel_tan      (DU/TU) tangential velocity
# |         Y(7) = mass         (kg) mass
# |
# |     -Used:
# |         [time,Y] = ode45(@MGALT_IN_FSM_2D_EOM,tspan,Y0,options,...
# |                     AU,TU,mew_sc)
# |
# |     -References
# |         B. A. Conway, Spacecraft Trajectory Optimization. Cambridge University Press, 2010.
# |
# |-----------------------------------------------------------------------

import math
import numpy as np
from Other_Functions.getSCThrust import getSCThrust
from Other_Functions.getSCmdot import getSCmdot

def MGALT_IN_FBSM_2D_EOM(t,Y,const,opts,mew_sc):

    ## EOM's
    
    # Unpack
    lamda1 = Y[0]          # (unitless) costate variable 1
    if lamda1 == 0:
        lamda1 = 0.000000000000001
    elif lamda1 > 1:
        lamda1 = 1
    elif lamda1 < -1:
        lamda1 = -1
    lamda2 = Y[1]          # (unitless) costate variable 2
    if lamda2 == 0:
        lamda2 = 0.000000000000001
    elif lamda2 > 1:
        lamda2 = 1
    elif lamda2 < -1:
        lamda2 = -1
    lamda3 = Y[2]          # (unitless) costate variable 3
    if lamda3 == 0:
        lamda3 = 0.000000000000001
    elif lamda3 > 1:
        lamda3 = 1
    elif lamda3 < -1:
        lamda3 = -1
    pos_rad = Y[3]         # (DU) radial position
    if pos_rad <= 0.000000001:
        pos_rad = 0.000000001
    # pos_ang = Y[4]     	   # (deg) angular position
    vel_rad = Y[5]      	# (DU/TU) radial velocity
    vel_tan = Y[6]       	# (DU/TU) tangential velocity
    mass = Y[7]          	# (kg) mass
    
    
    # Interpret if SC is thrusting and what the thrust is
    T = getSCThrust(const,opts,pos_rad,0)   # (kg*DU/TU^2)
     
    
    # Calculate Thrust Angle, Phi, With Costate Variables
    cosphi = -lamda2/math.sqrt((lamda1)**2+(lamda2)**2)   # (unitless) cos(phi) where phi is the thrust angle
    sinphi = -lamda1/math.sqrt((lamda1)**2+(lamda2)**2)   # (unitless) sin(phi) where phi is the thrust angle
    if np.isnan(cosphi) or np.isnan(sinphi):
        print("Lambda 1: %.15f" % lamda1)
        print("Lambda 2: %.15f" % lamda2)
        print('The damn thing broke again. Check these values ^^ to be positive. They are probably 0 or somehow negative.')
        1/0

    
    # Equations of Motion
    d_lamda1 = -lamda3 + vel_tan*lamda2/pos_rad
    d_lamda2 = (-2*vel_tan*lamda1+vel_rad*lamda2)/pos_rad
    d_lamda3 = (vel_tan**2*lamda1-vel_rad*vel_tan*lamda2)/pos_rad**2 - 2*mew_sc*lamda1/pos_rad**3
    d_pos_rad = vel_rad
    d_pos_ang = vel_tan/pos_rad*180/3.14159
    if np.isnan(d_pos_ang):
        print("Broke")
        1/0
    d_vel_rad = vel_tan**2/pos_rad - mew_sc/pos_rad**2 + T*sinphi/mass
    d_vel_tan = -vel_rad*vel_tan/pos_rad + T*cosphi/mass
    d_mass = getSCmdot(const,opts,T)       # kg/TU
    
    # Output
    dY = [d_lamda1,d_lamda2,d_lamda3,d_pos_rad,d_pos_ang,d_vel_rad,d_vel_tan,d_mass]

    return dY
