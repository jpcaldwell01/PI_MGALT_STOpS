# FORM: [R,V] = MGALT_stateBodies(JD_target,BOD,CONST,OPT,VAR,array_bodies)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Function to extract planetary location and velocity for certain 
# |     Julan Dates from large variable arrays
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -JD_target        	(1,1)       [float]         [JD]
# |         The desired JD to get planetary R and V at
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
# |     -array_bodies       (3,1)       [int]        [unitless]
# |         An array containing the index numbers for the bodies
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -R                  (6,1)       [float]         [AU]
# |         Radius vector [planet1 planet2]
# |     -V                  (6,1)       [float]         [AU/TU]
# |         Velocity vector [planet1 planet2]
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np
from scipy.integrate import odeint
from Other_Functions.orbit3D import orbit3D

def MGALT_stateBodies(JD_target,bod,const,opts,var,array_bodies):

    ## Get Info
    
    # Find which indicies the 0 value is located at
    index = np.argmin([abs(value) for value in bod['bodies_JD']-JD_target])
    is_backwards = (bod['bodies_JD'][index]-JD_target) >= 0
    JD_delta = bod['bodies_JD'][index]-JD_target          # How many JDs ahead of the target JD the closest segment is 
    
    # Check if JD_delta*86400 > -1 and < 1, due to tspan being indexed by 1
    if ((JD_delta*86400 > -1) and (JD_delta*86400 < 1)):
        
        R = [bod['bodies_R'][val][index] for val in array_bodies]
        V = [bod['bodies_V'][val][index] for val in array_bodies]
    
    # Check error for first or last index under/overflow
    if ((index == 1) and (JD_delta >= 1)) or ((index == len(bod['bodies_JD'])) and (JD_delta <= 1)):
        
        raise Exception('Error with JD parsing! The desired JD date is out of range of the valid JD dates stored in "BOD.bodies_JD". ODE45 will make solution inaccurate...exiting program')
    
    # If is_backwards, means backwards propogate to get to R and V values
    if is_backwards:
        tspan = list(range(0,round(-JD_delta*86400)-1,-1))
    else:
        tspan = list(range(0,round(-JD_delta*86400)+1,1))

    # Get R and V data
    R_init = [bod['bodies_R'][val][index] for val in array_bodies]
    V_init = [bod['bodies_V'][val][index] for val in array_bodies]
    state = R_init+V_init
    # Use solve_IVP to get the planet position at the desired JD
    # data = ode45(@orbit3D,tspan,[R_init;V_init],OPT.ode,const['Sun_mu']')[1]
    data = odeint(orbit3D,state,tspan,args=(const['Sun_mu'][0],),rtol=1e-8)
    
    # Final Results
    R = data[-1,0:3]
    V = data[-1,3:6]
    
    return R,V