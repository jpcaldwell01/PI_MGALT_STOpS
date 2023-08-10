#  FORM: [] = plotPareto(fs_trans_GA,fs_trans_DE,fs_trans_PSO,fs_trans_MBH,num_mig,opts)
# |
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Function for plotting the by-leg costs of a 2-transfer trajectory 
# |      over the course of optimization
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -fs_trans_GA    (opts['GA']['Npop']*opts['GA']['Ngen'],# of transfers)   [Array of Float]   [cost]
# |         An array carrying all individual transfer costs for the Genetic Algorithm Island. 
# |         Resets after each migration.
# |     -fs_trans_DE    (opts['DE']['Npop']*opts['DE']['Ngen'],# of transfers)   [Array of Float]   [cost]
# |         An array carrying all individual transfer costs for the Differential Evolution Island. 
# |         Resets after each migration.
# |     -fs_trans_PSO    (opts['PSO']['Npop']*opts['PSO']['tspan'],# of transfers)   [Array of Float]   [cost]
# |         An array carrying all individual transfer costs for the Particle Swarm Island. 
# |         Resets after each migration.
# |     -fs_trans_MBH    (Nkeep,# of transfers)   [Array of Float]   [cost]
# |         An array carrying all individual transfer costs for the Monotonic Basin Island.
# |         Resets after each migration.
# |     -num_mig             (1,1)       [int]         [unitless] 
# |         Number of migrations completed
# |     -opts                (1,1)       [dict]        [unitless]
# |         A struct containing constants user options. Contains the save 
# |         folder, ToF values, and more structs containing informaiton 
# |         for the island model, cost parameters, weighting parameters, 
# |         and all of the islands used in the optimization process
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------


import matplotlib.pyplot as plt

def plotPareto(fs_trans_GA,fs_trans_DE,fs_trans_PSO,fs_trans_MBH,num_mig,opts):
    fig = plt.figure()
    ax = fig.add_subplot()
    if 'GA' in opts['island']['isl_list']:
        ax.loglog(fs_trans_GA[:,0],fs_trans_GA[:,1],'rx',label='GA')
    if 'DE' in opts['island']['isl_list']:
        ax.loglog(fs_trans_DE[:,0],fs_trans_DE[:,1],'bx',label='DE')
    if 'PSO' in opts['island']['isl_list']:
        ax.loglog(fs_trans_PSO[:,0],fs_trans_PSO[:,1],'cx',label='PSO')
    if 'MBH' in opts['island']['isl_list']:
        ax.loglog(fs_trans_MBH[:,0],fs_trans_MBH[:,1],'gx',label='MBH')
    plt.title('Migration %i Island Costs' % num_mig)
    plt.xlabel('Cost E->M')
    plt.ylabel('Cost M->J')
    return