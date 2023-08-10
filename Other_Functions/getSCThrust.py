# FORM: [T] = getSCThrust(CONST,OPT,pos_rad,mem_thrust)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Function to get the current thrust value of the spacecraft. Used
# |     in the EOM ODE's
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
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
# |     -pos_rad            (1,1)       [float]         [AU]
# |         The radial position of the S/C with respect to the sun
# |     -mem_thrust         (1,n)       [boolean/float] [N]
# |         Optional, for the direct method to let STOpS know if there is
# |         thrust during that time segment (0 or 1), or if using variable 
# |         thrust, the amount of thrust
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -T                  (1,1)       [float]         [N]
# |         Spacecraft Thrust
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

def getSCThrust(const,opts,pos_rad,mem_thrust):

    ## Get SC Thrust
    
    if opts['thrust']['thrust_method'] == 'constant':       # Constant Thrust
        
        if opts['solver'] in ['LT_DIR_FSM_2D','MGALT_DIR_FBSM_2D']:
            
        # Check if the direct method is thrusting or not
            
            if mem_thrust:
                T = opts['thrust']['thrust']	# N
            else:
                T = 0                           # N
            
            
        elif opts['solver'] in ['LT_IN_FSM_2D','MGALT_IN_FBSM_2D']:
            
            T = opts['thrust']['thrust']      # N

        else:
            
            raise Exception('Invalid Solver')
        
    elif opts['thrust']['thrust_method'] == 'variable':
        
        if opts['solver'] in ['LT_DIR_FSM_2D','MGALT_DIR_FBSM_2D']:
                
            T = mem_thrust     # N
                
        elif opts['solver'] in ['LT_IN_FSM_2D','MGALT_IN_FBSM_2D']:
                
            raise Exception('The INDIRECT Method does not support varaible thrust.')
                
        else:
                
            raise Exception('Invalid thrust method.')
        
    elif opts['thrust']['thrust_method'] == 'equation':
        
        raise Exception('Equation Thrust Method unfinished')
        '''
        switch OPT.solver
            
            case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
                
                T = feval(OPT.thrust.thrust,OPT,pos_rad)       # N
                
            case {'LT_DIR_FSM_2D','MGALT_DIR_FBSM_2D'}
                
                if mem_thrust
                    T = feval(OPT.thrust.thrust,OPT,pos_rad)	# N
                else
                    T = 0                                      # N
                end
                
            otherwise
                
                errorPathDisplay()
                errorSolver()
                return
                
        end
        '''
        
    else:
        
        raise Exception('Incorrect thrust profile selected.')
    
    T = T/1000         # (kg*km/s^2)
    T = T*(const['TU'][0]*86400)**2/const['AU'][0]      # (kg*DU/TU^2)
    
    return T
