import PI_MGALT_STOpS
import PI_MGALT_STOpS_Demo

#%% Run Mission

# Run Demo
# [bod,opts,var,isBroken,best_sol] = PI_MGALT_STOpS_Demo.run(miss = 'demo_short') # 'demo_short' or 'demo_long'
                    # Do not expect the short demo to always produce good results.

# Run Custom Mission
[bod,opts,var,isBroken,best_sol] = PI_MGALT_STOpS.run(miss = 'MGALT_STOpS_T2_IN_SHORT')

#%% Post Processing
    

