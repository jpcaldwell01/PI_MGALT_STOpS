# FORM: [vel_sc_exit,rp] = MGALT_conditionsGravityAssist(BOD,CONST,...
#       pos_body,vel_body,vel_sc_enter,coe,index)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function gets the spacecraft velocity after a gravity assist,
# |     as well as the radius of perigee
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
# |     -pos_body           (3,1)    	[float]         [AU]
# |         Heliocentric position of the planet
# |     -vel_body           (3,1)       [float]         [AU/TU]
# |         Heliocentric velocity of the planet
# |     -vel_sc_enter      	(3,1)    	[float]         [AU/TU]
# |         Heliocentric velocity of the spacecraft before gravity assist
# |     -coe                (1,1)       [float]         [unitless]
# |         Coefficient to determine rp value
# |     -index            	(1,1)       [int]           [unitless]
# |         Used to determine current planet
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -vel_sc_exit      	(3,1)       [float]         [DU/TU]
# |         Heliocentric radial velocity of spacecraft after flyby
# |     -rp                 (1,1)       [float]         [km]
# |         Radius of perigee for spacecraft flyby
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np
import math

def MGALT_conditionsGravityAssist(bod,const,pos_body,vel_body,vel_sc_enter,coe,index):

    ## Initial
    
    # Some basic info
    pos_center = [0,0,0]
    mue_sun = const['Sun_mu'][0]                                 # km^3/s^2
    mue = const[bod['bodies'][index]+'_mu'][0]       	# km^3/s^2
    SOI = const[bod['bodies'][index]+'_SOI'][0]/2       # km *radius
    min_flyby = const[bod['bodies'][index]+"_rp"][0]	# km *radius
    
    
    ## Solve for the Exit Velocity Vectors

    # Solve for the unit vectors, Curtis: 8.75
    u_v = vel_body/np.linalg.norm(vel_body)              # planet's velocity vector
    planet_to_sun = pos_center - pos_body
    u_s = planet_to_sun/np.linalg.norm(planet_to_sun)	# planet pointing to the sun
    
    
    # Get the scalar components of V1_vec, Curtis: 8.76
    # Angle between V1_vec and V_vec, Curtis: 8.76
    cos_alpha = np.dot(vel_sc_enter,vel_body)/(np.linalg.norm(vel_sc_enter)*np.linalg.norm(vel_body))
    try:
        alpha_1 = math.acos(cos_alpha)                  # rad
    except:
        1/0
    
    
    # Calculate V, Curtis: 8.79
    V = math.sqrt(mue_sun/np.linalg.norm(pos_body))
    
    
    # Calculate V_inf1 vec, Curtis: 8.72
    V_inf1_vec = vel_sc_enter-vel_body
    V_inf1_mag = np.linalg.norm(V_inf1_vec)
    
    
    # Get the scalar components of V_inf1_vec, Curtis: 8.81
    V_inf1_V = np.linalg.norm(vel_sc_enter)*math.cos(alpha_1) - V
    V_inf1_S = np.linalg.norm(vel_sc_enter)*math.sin(alpha_1)
    
    
    # Angle between V_inf1_vec and V_vec, Curtis: 8.84
    phi_1 = math.atan2(V_inf1_S,V_inf1_V)   # rad
    if phi_1 < 0:    # Quadrant ambiguity
        phi_1 = (2*math.pi)+phi_1           # rad
    
    # Calculate the turning angle, delta, Curtis: 8.54
    # Get v_inf_mag, Curtis: 8.50 and 8.51
    if np.linalg.norm(vel_sc_enter) <= np.linalg.norm(vel_body):
        v_inf_mag = np.linalg.norm(vel_body)-np.linalg.norm(vel_sc_enter)
    else:
        v_inf_mag = np.linalg.norm(vel_sc_enter)-np.linalg.norm(vel_body)

    rp = coe*(SOI-min_flyby)+min_flyby  # Get radius perigee from un-normed data 
    
    delta = 2*math.asin(1/(1+((rp*(v_inf_mag**2))/mue)))     # rad
    
    
    # Angle between V_inf2_vec and V_vec, Curtis: 8.85
    phi_2 = phi_1+delta	# rad
    
    
    # Calculate V_inf2_vec, Curtis: 8.86
    V_inf2_vec = V_inf1_mag*math.cos(phi_2)*u_v + V_inf1_mag*math.sin(phi_2)*u_s
    
    
    # V2 vec, Curtis: 8.87
    vel_sc_exit = vel_body+V_inf2_vec

    return vel_sc_exit,rp
