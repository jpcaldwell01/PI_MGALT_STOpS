# FORM: [info] = MGALT_varLimits(BOD,OPT)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Function designed to generate variable limits and binary check 
# |     array for low-thrust trajectory variable strings
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -BOD                (1,1)       [struct]        [unitless]
# |         A struct containing information pertaining to the planetary
# |         bodies. Contains list of bodies, launch windows and ToF, and 
# |         planetary R/V/JD vectors. This struct has dynamic fields and 
# |         will adapt to contain only the necesary information
# |     -OPT                (1,1)       [struct]        [unitless]
# |         A struct containing constants user options. Contains the save 
# |         folder, ToF values, and more structs containing informaiton 
# |         for the island model, cost parameters, weighting parameters, 
# |         and all of the islands used in the optimization process
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -info               (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np
from julian import to_jd

def MGALT_varLimits(bod,opts):
    
    #%% Setup
    
    # opt_thrust = OPT.thrust
    bodies = bod['bodies']
    
    # Desired transfer times for each segment
    Htime = opts['thrust']['tt_end']
    transfers = len(bodies)-1
    
    
    
    ### Method Selection
    if opts['solver'] == 'LT_DIR_FSM_2D':	# Segmented Method

        # Preallocate: Lows/Highs/Binary
        low = np.zeros((transfers,opts['thrust']['Nseg']*2+2))
        high = np.zeros((transfers,opts['thrust']['Nseg']*2+2))
        bins = np.zeros((transfers,opts['thrust']['Nseg']*2+2))
        
        
        # Departure body
        low[0,0] = to_jd(bod['window1'])      # earliest departure time (JD)
        high[0,0] = to_jd(bod['window2'])     # latest departure time (JD)
        
        if opts['thrust']['thrust_method'] in ['constant','equation']:
            
            # This thrust method assumes that "first Tswitch assumed 1"
            # from Sheehan's code base. It was kept the same way while
            # updating to MGALT, but only for the departure planet.
            # It was assumed only for the departure planet because
            # subsequent flybys may need the thrust to be off after
            # happening
            
            # Departure body
            low[0,1] = 1   # Tswitch low (on)
            high[0,1] = 1  # Tswitch high (on)
            bins[0,1] = 0   # non-binary variable, to force the thrust to always be on for departure

            # Segments between departure and next planet
            for i1 in range(2,opts['thrust']['Nseg']*2+1,2):     
                low[0,i1] = 0      # phi low
                high[0,i1] = 360   # phi high
                bins[0,i1] = 0      # non-binary variable
                
                low[0,i1+1] = 0  	# Tswitch low (off)
                high[0,i1+1] = 1 	# Tswitch high (on)
                bins[0,i1+1] = 1   	# binary variable  
            
            low[0,-1] = Htime[0] - opts['thrust']['time'][0][0]	# shortest transfer time (days)
            high[0,-1] = Htime[0] + opts['thrust']['time'][0][1]	# longest transfer time (days)
            bins[0,-1] = 0                                 # non-binary variable
            
            # Other transfer segment departures
            for i2 in range(1,transfers):
                
                low[i2] = low[0]
                high[i2] = high[0]
                bins[i2] = bins[0]
                
                low[i2,0] = low[i2-1,0]+low[i2-1,-1]              # earliest departure time
                high[i2,0] = high[i2-1,0]+high[i2-1,-1]           # latest departure time
                low[i2,-1] = Htime[i2] - opts['thrust']['time'][i2][0]	# shortest transfer time (days)
                high[i2,-1] = Htime[i2] + opts['thrust']['time'][i2][1]	# longest transfer time (days)
                
            
            # Make the lows for every transfer be 0 so it's not the
            # assumption of having a powered thrust like from the
            # departure body
            low[1:,1] = 0
            bins[1:,1] = 1   # For having preliminary segments be powered/unpowered after transfer
         
        elif  opts['thrust']['thrust_method'] == 'variable':
            
            # Departure body
            for i1 in range(0,opts['thrust']['Nseg']*2,2):
                low[0,i1+1] = opts['thrust']['thrust'][0]     # T low
                high[0,i1+1] = opts['thrust']['thrust'][1]    # T high
                
                low[0,i1+2] = 0                        # phi low
                high[0,i1+2] = 360                     # phi high
            
            low[0,-1] = Htime[0] - opts['thrust']['time'][0][0]   # shortest transfer time (days)
            high[0,-1] = Htime[0] + opts['thrust']['time'][0][1]  # longest transfer time (days)
            bins[0,-1] = 0                                 # non-binary variable
            
            #Other transfer segment departures
            for i2 in range(1,transfers):
                
                low[i2] = low[0]
                high[i2] = high[0]
                
                low[i2,0] = low[i2-1,0]+ low[i2-1,-1]      # earliest departure time
                high[i2,0] = high[i2-1,0] + high[i2-1,-1]	# latest departure time
                low[i2,-1] = Htime[i2] - opts['thrust']['time'][0,0]    # shortest transfer time (days)
                high[i2,-1] = Htime[i2] + opts['thrust']['time'][0,1]	# longest transfer time (days)
            
        else:
            
            raise Exception('Invalid Thrust Method')
            return
      
    elif opts['solver'] == 'LT_IN_FSM_2D':    	# Costate Method

        # Preallocate: Lows/Highs/Binary
        low = np.zeros((transfers,5))
        high = np.zeros((transfers,5))
        bins = np.zeros((transfers,5))

        # Departure body
        low[0,0] = to_jd(bod['window1'])          # earliest departure time (JD)
        high[0,0] = to_jd(bod['window2'])         # latest departure time (JD)
        low[0,1:4] = -1                                # minimum of -1 for lamda variables
        high[0,1:4] = 1                                # maximium of 1 for lamda variables
        low[0,4] = Htime[0] - opts['thrust']['time'][0][0]     # shortest transfer time (days)
        high[0,4] = Htime[0] + opts['thrust']['time'][0][1]	# longest transfer time (days) 
        
        #Other transfer segment departures
        for i1 in range(1,transfers):
            low[i1,0] = low[i1-1,0] + low[i1-1,4]            # earliest departure time (JD)
            high[i1,0] = high[i1-1,0] + high[i1-1,4]         # latest departure time (JD)
            low[i1,1:4] = -1                               # minimum of -1 for lamda variables
            high[i1,1:4] = 1                               # maximium of 1 for lamda variables
            low[i1,4] = Htime(i1) - opts['thrust']['time'][0][0]  # shortest transfer time (days)
            high[i1,4] = Htime(i1) + opts['thrust']['time'][0][1] # longest transfer time (days)     
            
    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':	# Segmented Method
        
        # Make sure segments are even
        if opts['thrust']['Nseg'] % 2:
            
            raise Exception('The DIRECT FBSM only supports an even number of segments')
            return
        
        # Preallocate: Lows/Highs/Binary
        low_dep = np.zeros((1,opts['thrust']['Nseg']+1))
        high_dep = np.zeros((1,opts['thrust']['Nseg']+1))
        bin_dep = np.zeros((1,opts['thrust']['Nseg']+1))
        
        low_ari = np.zeros((1,opts['thrust']['Nseg']+1))
        high_ari = np.zeros((1,opts['thrust']['Nseg']+1))
        bin_ari = np.zeros((1,opts['thrust']['Nseg']+1))

        low_trans = np.zeros((transfers-1,2*opts['thrust']['Nseg']+5))
        high_trans = np.zeros((transfers-1,2*opts['thrust']['Nseg']+5))
        bin_trans = np.zeros((transfers-1,2*opts['thrust']['Nseg']+5))

        reshaped_low_trans = np.array([])
        reshaped_high_trans = np.array([])
        reshaped_bin_trans = np.array([])

        # Departure body
        low_dep[0,0] = to_jd(bod['window1'])    	# earliest departure time (JD)
        high_dep[0,0] = to_jd(bod['window2'])     # latest departure time (JD)
        bin_dep[0,0] = 0                               # non-binary variable
        
        
        # Thrust methods        
        
        if opts['thrust']['thrust_method'] in ['constant','equation']:
            '''
            Add constraint pruning
            # This thrust method assumes that "first Tswitch assumed 1"
            # from Sheehan's code base. It was kept the same way while
            # updating to MGALT, but only for the departure planet.
            # It was assumed only for the departure planet because
            # subsequent flybys may need the thrust to be off after
            # happening
            
            # Departure body
            low_dep[0,1] = 1   # Tswitch low (on)
            high_dep[0,1] = 1  # Tswitch high (on)
            bin_dep[0,1] = 0   # non-binary variable, to force the thrust to always be on for departure

            # Departure segments
            for i1 in range(3,opts['thrust']['Nseg'],2):
                low_dep[0,i1-1] = 0      # phi low
                high_dep[0,i1-1] = 360   # phi high
                bin_dep[0,i1-1] = 0
                
                low_dep[0,i1] = 0  	# Tswitch low (off)
                high_dep[0,i1] = 1 	# Tswitch high (on)
                bin_dep[0,i1] = 1   	# binary variable  
            
            low_dep[0,opts['thrust']['Nseg']] = 0
            high_dep[0,opts['thrust']['Nseg']] = 360
            bin_dep[0,opts['thrust']['Nseg']] = 0
            
            # Arrival segments
            for i2 in range(1,opts['thrust']['Nseg'],2):
                low_ari[0,i2-1] = 0      # Tswitch low (off)
                high_ari[0,i2-1] = 1     # Tswitch high (on)
                bin_ari[0,i2-1] = 1      # binary variable 
                
                low_ari[0,i2] = 0    # phi low
                high_ari[0,i2] = 360 # phi high
                bin_ari[0,i2] = 0
                
            low_ari[0,opts['thrust']['Nseg']] = Htime[-1] - opts['thrust']['time'][-1][0]     # earliest departure time
            high_ari[0,opts['thrust']['Nseg']] = Htime[-1] + opts['thrust']['time'][-1][1] 	# latest departure time
            bin_ari[0,opts['thrust']['Nseg']] = 0                                       # non-binary variable
            
            
            # For transfer segments
            if transfers > 1:
                
                # All transfers
                for i3 in range(transfers-1):
                    
                    # Prior to Flyby
                    for i4 in range(0,opts['thrust']['Nseg'],2):
                        low_trans[i3,i4] = 0       # Tswitch low (off)
                        high_trans[i3,i4] = 1      # Tswitch high (on)
                        bin_trans[i3,i4] = 1       # binary variable
                        
                        low_trans[i3,i4+1] = 0    	# phi low
                        high_trans[i3,i4+1] = 360  # phi high  
                        bin_trans[i3,i4+1] = 0
                
                    # Flyby Variables
                    low_trans[i3,opts['thrust']['Nseg']] = 0                        # minimum of 0 for the flyby coefficient for rp
                    high_trans[i3,opts['thrust']['Nseg']] = 1                       # maximum of 1 for the flyby coefficient for rp
                    low_trans[i3,opts['thrust']['Nseg']+1:opts['thrust']['Nseg']+3] = -1     # Englander up to unit vector control min
                    high_trans[i3,opts['thrust']['Nseg']+1:opts['thrust']['Nseg']+3] = 1     # Englander up to unit vector control max
                    low_trans[i3,opts['thrust']['Nseg']+3] = Htime[i3] - opts['thrust']['time'][i3][0]                        # shortest transfer time (days)
                    high_trans[i3,opts['thrust']['Nseg']+3] = Htime[i3] + opts['thrust']['time'][i3][1]                        # longest transfer time (days) 
                    low_trans[i3,opts['thrust']['Nseg']+4] = low_dep[0,0] + sum(low_trans.transpose()[opts['thrust']['Nseg']+3][:])   # earliest departure time (JD)
                    high_trans[i3,opts['thrust']['Nseg']+4] = high_dep[0,0] + sum(high_trans.transpose()[opts['thrust']['Nseg']+3][:])	# latest departure time (JD)
                    bin_trans[i3,opts['thrust']['Nseg']:opts['thrust']['Nseg']+5] = 0
                    
                    # After Flyby
                    for i5 in range(opts['thrust']['Nseg']+5,2*opts['thrust']['Nseg']+4,2):
                        low_trans[i3,i5] = 0       # Tswitch low (off)
                        high_trans[i3,i5] = 1      # Tswitch high (on)
                        bin_trans[i3,i5] = 1       # binary variable
                        
                        low_trans[i3,i5+1] = 0    	# phi low
                        high_trans[i3,i5+1] = 360	# phi high
                        bin_trans[i3,i5+1] = 0       # binary variable
                
                # Reshape the vectors          
                for i6 in range(transfers-1):
                    reshaped_low_trans = np.append(reshaped_low_trans, low_trans[i6])
                    reshaped_high_trans = np.append(reshaped_high_trans, high_trans[i6])
                    reshaped_bin_trans = np.append(reshaped_bin_trans, bin_trans[i6])

        '''
        elif opts['thrust']['thrust_method'] ==  'variable':
            
            # This thrust method assumes that "first Tswitch assumed 1"
            # from Sheehan's code base. It was kept the same way while
            # updating to MGALT, but only for the departure planet.
            # It was assumed only for the departure planet because
            # subsequent flybys may need the thrust to be off after
            # happening
            
            # Departure body
            low_dep[0,1] = opts['thrust']['thrust'][0]  	# T low
            high_dep[0,1] = opts['thrust']['thrust'][1]   # T high
            bin_dep = np.zeros((1,opts['thrust']['Nseg']+1))   # non-binary variable
            
            if opts['constraint_Pruning'] == 'n':
                # Departure segments
                for i1 in range(2,opts['thrust']['Nseg'],2):
                    
                    low_dep[0,i1] = 0      # phi low
                    high_dep[0,i1] = 360   # phi high
                    
                    low_dep[0,i1+1] = opts['thrust']['thrust'][0]     # T low
                    high_dep[0,i1+1] = opts['thrust']['thrust'][1]	# T high
                    
                
                low_dep[0,opts['thrust']['Nseg']] = 0
                high_dep[0,opts['thrust']['Nseg']] = 360
                
                # Arrival segments
                for i2 in range(0,opts['thrust']['Nseg'],2):
                    low_ari[0,i2] = opts['thrust']['thrust'][0]	# T low
                    high_ari[0,i2] = opts['thrust']['thrust'][1]	# T high
                    
                    low_ari[0,i2+1] = 0    # phi low
                    high_ari[0,i2+1] = 360 # phi high
            
            elif opts['constraint_Pruning'] == 'y':
                constraintRange = [45,75] # (Both Positive) Constrains the initial thruster angle to be within -/+ degrees of the planet's velocity
                if np.linalg.norm([bod['bodies_R'][0][0],bod['bodies_R'][1][0],bod['bodies_R'][2][0]]) < np.linalg.norm([bod['bodies_R'][3][0],bod['bodies_R'][4][0],bod['bodies_R'][5][0]]):
                    # Departure segments
                    for i1 in range(2,opts['thrust']['Nseg'],2):
                        
                        low_dep[0,i1] = -constraintRange[0]      # phi low
                        high_dep[0,i1] = constraintRange[1]   # phi high
                        
                        low_dep[0,i1+1] = opts['thrust']['thrust'][0]     # T low
                        high_dep[0,i1+1] = opts['thrust']['thrust'][1]	# T high
                        
                    
                    low_dep[0,opts['thrust']['Nseg']] = -constraintRange[0]
                    high_dep[0,opts['thrust']['Nseg']] = constraintRange[1]
                    
                    # Arrival segments
                    for i2 in range(0,opts['thrust']['Nseg'],2):
                        low_ari[0,i2] = opts['thrust']['thrust'][0]	# T low
                        high_ari[0,i2] = opts['thrust']['thrust'][1]	# T high
                        
                        low_ari[0,i2+1] = -constraintRange[0]    # phi low
                        high_ari[0,i2+1] = constraintRange[1] # phi high
                
                elif np.linalg.norm([bod['bodies_R'][0][0],bod['bodies_R'][1][0],bod['bodies_R'][2][0]]) >= np.linalg.norm([bod['bodies_R'][3][0],bod['bodies_R'][4][0],bod['bodies_R'][5][0]]):
                    # Departure segments
                    for i1 in range(2,opts['thrust']['Nseg'],2):
                        
                        low_dep[0,i1] = -constraintRange[0] + 180      # phi low
                        high_dep[0,i1] = constraintRange[1] + 180   # phi high
                        
                        low_dep[0,i1+1] = opts['thrust']['thrust'][0]     # T low
                        high_dep[0,i1+1] = opts['thrust']['thrust'][1]	# T high
                        
                    
                    low_dep[0,opts['thrust']['Nseg']] = -constraintRange[0] + 180
                    high_dep[0,opts['thrust']['Nseg']] = constraintRange[1] + 180
                    
                    # Arrival segments
                    for i2 in range(0,opts['thrust']['Nseg'],2):
                        low_ari[0,i2] = opts['thrust']['thrust'][0]	# T low
                        high_ari[0,i2] = opts['thrust']['thrust'][1]	# T high
                        
                        low_ari[0,i2+1] = -constraintRange[0] + 180    # phi low
                        high_ari[0,i2+1] = constraintRange[1] + 180 # phi high
                        
            low_ari[0,opts['thrust']['Nseg']] = Htime[-1] - opts['thrust']['time'][-1][0]    # earliest departure time
            high_ari[0,opts['thrust']['Nseg']] = Htime[-1] + opts['thrust']['time'][-1][1]	# latest departure time
            
            
            # For transfer segments
            if transfers > 1:
                if opts['constraint_Pruning'] == 'n':
                    
                    # All transfers
                    for i3 in range(transfers-1):
                        
                        # Prior to Flyby
                        for i4 in range(0,opts['thrust']['Nseg'],2):
    
                            low_trans[i3,i4] = opts['thrust']['thrust'][0]	# T low
                            high_trans[i3,i4] = opts['thrust']['thrust'][1]	# T high
                            
                            low_trans[i3,i4+1] = 0    	# phi low
                            high_trans[i3,i4+1] = 360  # phi high                
                        
                        # Flyby Variables
                        low_trans[i3,opts['thrust']['Nseg']] = 0                        # minimum of 0 for the flyby coefficient for rp
                        high_trans[i3,opts['thrust']['Nseg']] = 1                       # maximum of 1 for the flyby coefficient for rp
                        low_trans[i3,opts['thrust']['Nseg']+1:opts['thrust']['Nseg']+3] = -1     # Englander up to unit vector control min
                        high_trans[i3,opts['thrust']['Nseg']+1:opts['thrust']['Nseg']+3] = 1     # Englander up to unit vector control max
                        low_trans[i3,opts['thrust']['Nseg']+3] = Htime[i3] - opts['thrust']['time'][i3][0]                        # shortest transfer time (days)
                        high_trans[i3,opts['thrust']['Nseg']+3] = Htime[i3] + opts['thrust']['time'][i3][1]                       # longest transfer time (days) 
                        low_trans[i3,opts['thrust']['Nseg']+4] = low_dep[0,0] + sum(low_trans.transpose()[opts['thrust']['Nseg']+3])       # earliest departure time (JD)
                        high_trans[i3,opts['thrust']['Nseg']+4] = high_dep[0,0] + sum(high_trans.transpose()[opts['thrust']['Nseg']+3])	# latest departure time (JD)
                        
                        # After Flyby
                        for i5 in range(opts['thrust']['Nseg']+5,2*opts['thrust']['Nseg']+4,2):
                            low_trans[i3,i5] = opts['thrust']['thrust'][0]	# T low
                            high_trans[i3,i5] = opts['thrust']['thrust'][1]	# T high
                            
                            low_trans[i3,i5+1] = 0    	# phi low
                            high_trans[i3,i5+1] = 360	# phi high
                            
                    # Reshape the vectors          
                    for i6 in range(transfers-1):
                        reshaped_low_trans = np.append(reshaped_low_trans, low_trans[i6])
                        reshaped_high_trans = np.append(reshaped_high_trans, high_trans[i6])
                        reshaped_bin_trans = np.append(reshaped_bin_trans, bin_trans[i6])
                           
                elif opts['constraint_Pruning'] == 'y':
                    # All transfers
                    for i3 in range(transfers-1):
                        if np.linalg.norm([bod['bodies_R'][3*(i3)][0],bod['bodies_R'][3*(i3)+1][0],bod['bodies_R'][3*(i3)+2][0]]) >= np.linalg.norm([bod['bodies_R'][3*(i3+1)][0],bod['bodies_R'][3*(i3+1)+1][0],bod['bodies_R'][3*(i3+1)+2][0]]):

                            # Prior to Flyby
                            for i4 in range(0,opts['thrust']['Nseg'],2):
        
                                low_trans[i3,i4] = opts['thrust']['thrust'][0]	# T low
                                high_trans[i3,i4] = opts['thrust']['thrust'][1]	# T high
                                
                                low_trans[i3,i4+1] = -constraintRange[0] + 180    	# phi low
                                high_trans[i3,i4+1] = constraintRange[1] + 180  # phi high  
                        
                        elif np.linalg.norm([bod['bodies_R'][3*(i3)][0],bod['bodies_R'][3*(i3)+1][0],bod['bodies_R'][3*(i3)+2][0]]) < np.linalg.norm([bod['bodies_R'][3*(i3+1)][0],bod['bodies_R'][3*(i3+1)+1][0],bod['bodies_R'][3*(i3+1)+2][0]]):
                            # Prior to Flyby
                            for i4 in range(0,opts['thrust']['Nseg'],2):
        
                                low_trans[i3,i4] = opts['thrust']['thrust'][0]	# T low
                                high_trans[i3,i4] = opts['thrust']['thrust'][1]	# T high
                                
                                low_trans[i3,i4+1] = -constraintRange[0]   	# phi low
                                high_trans[i3,i4+1] = constraintRange[1]    # phi high   
                                
                        # Flyby Variables
                        low_trans[i3,opts['thrust']['Nseg']] = 0                        # minimum of 0 for the flyby coefficient for rp
                        high_trans[i3,opts['thrust']['Nseg']] = 1                       # maximum of 1 for the flyby coefficient for rp
                        low_trans[i3,opts['thrust']['Nseg']+1:opts['thrust']['Nseg']+3] = -1     # Englander up to unit vector control min
                        high_trans[i3,opts['thrust']['Nseg']+1:opts['thrust']['Nseg']+3] = 1     # Englander up to unit vector control max
                        low_trans[i3,opts['thrust']['Nseg']+3] = Htime[i3] - opts['thrust']['time'][i3][0]                        # shortest transfer time (days)
                        high_trans[i3,opts['thrust']['Nseg']+3] = Htime[i3] + opts['thrust']['time'][i3][1]                       # longest transfer time (days) 
                        low_trans[i3,opts['thrust']['Nseg']+4] = low_dep[0,0] + sum(low_trans.transpose()[opts['thrust']['Nseg']+3])       # earliest departure time (JD)
                        high_trans[i3,opts['thrust']['Nseg']+4] = high_dep[0,0] + sum(high_trans.transpose()[opts['thrust']['Nseg']+3])	# latest departure time (JD)
                        
                        if np.linalg.norm([bod['bodies_R'][3*(i3+1)][0],bod['bodies_R'][3*(i3+1)+1][0],bod['bodies_R'][3*(i3+1)+2][0]]) < np.linalg.norm([bod['bodies_R'][3*(i3+2)][0],bod['bodies_R'][3*(i3+2)+1][0],bod['bodies_R'][3*(i3+2)+2][0]]):
                            # After Flyby
                            for i5 in range(opts['thrust']['Nseg']+5,2*opts['thrust']['Nseg']+4,2):
                                low_trans[i3,i5] = opts['thrust']['thrust'][0]	# T low
                                high_trans[i3,i5] = opts['thrust']['thrust'][1]	# T high
                                
                                low_trans[i3,i5+1] = -constraintRange[0]   	# phi low
                                high_trans[i3,i5+1] = constraintRange[1]	# phi high
                                
                        elif np.linalg.norm([bod['bodies_R'][3*(i3+1)][0],bod['bodies_R'][3*(i3+1)+1][0],bod['bodies_R'][3*(i3+1)+2][0]]) >= np.linalg.norm([bod['bodies_R'][3*(i3+2)][0],bod['bodies_R'][3*(i3+2)+1][0],bod['bodies_R'][3*(i3+2)+2][0]]):
                            # After Flyby
                            for i5 in range(opts['thrust']['Nseg']+5,2*opts['thrust']['Nseg']+4,2):
                                low_trans[i3,i5] = opts['thrust']['thrust'][0]	# T low
                                high_trans[i3,i5] = opts['thrust']['thrust'][1]	# T high
                                
                                low_trans[i3,i5+1] = -constraintRange[0] + 180   	# phi low
                                high_trans[i3,i5+1] = constraintRange[1] + 180	# phi high
                            
                    # Reshape the vectors          
                    for i6 in range(transfers-1):
                        reshaped_low_trans = np.append(reshaped_low_trans, low_trans[i6])
                        reshaped_high_trans = np.append(reshaped_high_trans, high_trans[i6])
                        reshaped_bin_trans = np.append(reshaped_bin_trans, bin_trans[i6])
        else:
            
            raise Exception('Invalid Thrust Method')
        
    elif opts['solver'] == 'MGALT_IN_FBSM_2D':    	# Costate Method
        
        # Preallocate: Lows/Highs/Binary
        low_trans = np.zeros((transfers-1,11))
        high_trans = np.zeros((transfers-1,11))
        bin_trans = np.zeros((transfers-1,11))

        reshaped_low_trans = np.empty((1,0))
        reshaped_high_trans = np.empty((1,0))
        reshaped_bin_trans = np.empty((1,0))

        
        # Departure body
        low_dep = [[]*4]
        high_dep = [[]*4]
        bin_dep = []
        
        low_dep[0] = to_jd(bod['window1'])    	# earliest departure time (JD)
        high_dep[0] = to_jd(bod['window2'])     # latest departure time (JD)
        low_dep[1:3] = [-1]*3                         	# minimum of -1 for dep lamda variables
        high_dep[1:3] = [1]*3                           	# maximium of 1 for dep lamda variables
        bin_dep[0:3] = [0]*4                             # non-binary variable
        
        # Arrival Body
        low_ari = [[]*4]
        high_ari = [[]*4]
        bin_ari = []
        
        low_ari[0:3] = [-1]*4                                    # minimum of -1 for ari lamda variables
        high_ari[0:3] = [1]*4                                    # maximium of 1 for ari lamda variables
        low_ari[3] = Htime[-1] - opts['thrust']['time'][-1][0]     # shortest transfer time (days)
        high_ari[3] = Htime[-1] + opts['thrust']['time'][-1][1]    # longest transfer time (days)
        bin_ari[0:3] = [0]*4                                     # non-binary variable
        
        # If there are transfers
        if transfers > 1:
            
            # Transfer segments
            for i1 in range(transfers-1):
                bin_trans[i1,0:10] = [0]*10     # non-binary variable
                low_trans[i1,0:3] = [-1]*3     # minimum of -1 for forward lamda variables
                high_trans[i1,0:3] = [1]*3     # maximum of 1 for forward lamda variables
                low_trans[i1,3] = 0        # minimum of 0 for the flyby coefficient for rp
                high_trans[i1,3] = 1       # maximum of 1 for the flyby coefficient for rp
                low_trans[i1,4:6] = [-1]*2     # Englander up to unit vector control min
                high_trans[i1,4:6] = [1]*2     # Englander up to unit vector control max
                low_trans[i1,6] = Htime[i1] - opts['thrust']['time'][i1][0]    # shortest transfer time (days)
                high_trans[i1,6] = Htime[i1] + opts['thrust']['time'][i1][1]   # longest transfer time (days) 
                low_trans[i1,7] = low_dep[0]+ sum(low_trans.transpose()[6][:])         # earliest departure time (JD)
                high_trans[i1,7] = high_dep[0]+ sum(high_trans.transpose()[6][:])      # latest departure time (JD)
                low_trans[i1,8:11] = [-1]*3     # minimum of -1 for backward lamda variables
                high_trans[i1,8:11] = [1]*3     # maximum of 1 for backward lamda variables
            
            
            # Reshape the vectors          
            for i2 in range(transfers-1):
                reshaped_low_trans = np.append(reshaped_low_trans, low_trans[i2])
                reshaped_high_trans = np.append(reshaped_high_trans, high_trans[i2])
                reshaped_bin_trans = np.append(reshaped_bin_trans, bin_trans[i2])
    else:                               # Error
        
        print('Incorrect solver selected.')
        print('Valid options are:')
        print('LT_DIR_FSM_2D')
        print('LT_IN_FSM_2D')
        print('MGALT_DIR_FBSM_2D')
        print('MGALT_IN_FBSM_2D')
        raise Exception("Please choose one of the correct solvers.")
        return
    
    ## Reshape the size
    
    reshaped_low = np.array([])
    reshaped_high = np.array([])
    reshaped_bin = np.array([])

    if opts['solver'] in ['LT_DIR_FSM_2D', 'LT_IN_FSM_2D']:
        
        reshaped_low = low[0]
        reshaped_high = high[0]
        reshaped_bin = bins[0]
        
        for i1 in range(1,transfers):
            reshaped_low = np.append(reshaped_low, low[i1])
            reshaped_high = np.append(reshaped_high, high[i1])
            reshaped_bin = np.append(reshaped_bin, bins[i1])

    elif  opts['solver'] in ['MGALT_DIR_FBSM_2D', 'MGALT_IN_FBSM_2D']:
        
        reshaped_low = np.append(np.append(np.append(reshaped_low,low_dep), reshaped_low_trans), low_ari)
        reshaped_high = np.append(np.append(np.append(reshaped_high,high_dep), reshaped_high_trans), high_ari)
        reshaped_bin = np.append(np.append(np.append(reshaped_bin,bin_dep), reshaped_bin_trans), bin_ari)
    
    info = {}    
    info['low'] = reshaped_low
    info['high'] = reshaped_high
    info['lowC'] = 1*info['low']
    info['highC'] = 1*info['high']
    info['bin'] = reshaped_bin
    info['transfers'] = transfers
    
    info['bestLegCosts'] = {}
    info['bestLegMembers'] = {}
    info['ps'] = .05
    info['bestLegsKept'] = 100
    for it in range(transfers):
        info['bestLegCosts']['leg%i' %(it)] = np.ndarray((0,transfers))
        info['bestLegMembers']['leg%i' %(it)] = np.ndarray((0,len(info['low'])))
    
    
    return info
