# FORM: [times,thrust,phi] = DIR_sepVariables(solver,member,index,seg)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Function designed to parse time, thrust, and thrust pointing 
# |     angle arrays from the MGALT_DIR_FSM_2D method
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -solver             (1,n)       [string]        [unitless]
# |         The solver (cost function) handle which is used to determine 
# |         the correct solver to use. This string is used in an ode45 
# |     -member             (1,Nvar)    [float]         [unitless]
# |         A single member of a population input into the function
# |     -index              (1,1)       [int]       	[unitless]
# |         The index value of a for loop used as an input to other 
# |         functions. The index value is useful in etracting information 
# |         like the planet string, R/V arrays, etc... without having to 
# |         pass in extra/unuseful information
# |     -seg             	(1,1)       [int]           [unitless]
# |         The segments for the direct method 
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -times              (1,2)       [float]         [JD][Day]
# |         The JD and ToF parsed from the member array
# |     -thrust             (1,Nseg)	[bool][float] 	[unitless][N]
# |         Binary indicator if the s/c is thrusting or not
# |         The thrust for each angle
# |     -phi                (1,Nseg) 	[float]         [deg]
# |         Thrust poointing angle for each segment
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------
import numpy as np

def DIR_sepVariables(solver,member,index,seg):

    ## Sep Vars
    times = np.zeros((1,2))
    seg = int(seg)
    
    if solver == 'LT_DIR_FSM_2D':
    
        times[0,0] = member[0]
        times[0,1] = member[-1]

        count = 0
        thrust = np.zeros((1,len(list(range(1,len(member)-1,2)))))
        phi = np.zeros((1,len(list(range(1,len(member)-1,2)))))

        for i1 in range(1,len(member)-1,2):
            thrust[0,count] = member[i1]   # N
            phi[0,count] = member[i1+1]    # deg
            count = count + 1
            
    elif solver == 'MGALT_DIR_FBSM_2D':
            
        # Only call this in the for i1 = 2:transfers-1 portion, as the
        # start and end locations are different
        times[0,0] = member[index*((2*seg)+5)]
        times[0,1] = member[(index+1)*((2*seg)+5)-1]
        data = member[   ((index*((2*seg)+5))+1) : ((index+1)*((2*seg)+5)-4)   ]
        thrust = np.array([data[it] for it in range(0,len(data),2)])
        thrust.shape = (1,len(list(range(0,len(data),2))))
        phi = np.array([data[it] for it in range(1,len(data),2)])
        phi.shape = (1,len(list(range(0,len(data),2))))
        
    else:
        
        raise Exception("Invalid Solver Selection")
    
    return times,thrust,phi
