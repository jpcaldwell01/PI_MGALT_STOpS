# FORM: [popn] = ISL_modReplacement(OPT_island,selected,num_mig,num_isl,popn)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Replacement function for solution sharing between islands. 
# |     This function selects the desired solutions for an island to 
# |     accept from the list of solutions that were shared by other 
# |     islands
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -OPT_island         (1,1)       [struct]        [unitless]
# |         Island Model option parameters
# |     -selected          	(Nmig,Nisl) [struct]        [unitless]
# |         Shared solutions from each island for each migration
# |     -num_mig            (1,1)       [int]           [unitless]
# |         The current migration number
# |     -num_isl            (1,1)       [int]           [unitless]
# |         The current island number
# |     -popn           (Nshared,Nvar)	[float]         [unitless]
# |         A randomely generated array of members from the initalization 
# |         of an island. Shared solutions will replace some or all of the 
# |         solutions in this array
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -popn           (Nshared,Nvar)	[float]         [unitless]
# |         The inital array or random members with the accepted 
# |         shared solutions added
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np

def ISL_modReplacement(opt_island,selected,num_mig,num_isl,popn):
    
    ## Setup
    
    [Npop,Nvar] = np.shape(popn)
    con_mat = opt_island['isl_conn'] # (giver, receiver) - connection matrix
    
    ## Check To See If Soultions Were Shared
    
    if len(selected[num_mig][num_isl]['f']) != 0:
       
    	# Find How Many Islands Shared
        n_share_isl = len(selected[num_mig])
       
        n_share_soln = np.zeros((n_share_isl,1))
        # Find How Many Solutions Each Island Shared
        for s_i in range(n_share_isl):
            n_share_soln[s_i] = np.shape(selected[num_mig][s_i]['sol'])[0]
        
        # Preallocate
        shared_soln = np.zeros((sum(con_mat[num_mig])*opt_island['sel_opt'][num_isl],Nvar))
        shared_f = np.zeros((sum(con_mat[num_mig])*opt_island['sel_opt'][num_isl],1))
       
        # Create Array of Shared Solutions - Only Include From Compatable Givers
        rep_sol = 0
        for s_i in range(n_share_isl):
            for s_s in range(int(n_share_soln[s_i])):
                # Verify Connection with Connection Matrix
                if con_mat[s_i][num_isl]:
                    shared_soln[rep_sol,:] = selected[num_mig][s_i]['sol'][s_s,:]
                    shared_f[rep_sol] = selected[num_mig][s_i]['f'][s_s]
                    rep_sol = rep_sol + 1
       
        # Sort The Shared Solutions
        sorted_f = np.sort(shared_f,axis=0)
        I = np.argsort(shared_f,axis=0)
        sorted_shared_soln = shared_soln[I][:,0,:]
        
        # Switch-Case based off of the sharing model
        if opt_island['rep_pol'][num_isl] == 'all':                # All
            '''  
            for s_s = 1:length(sorted_f)
                if s_s > Npop
                    break
                end
                popn(s_s,:) = sorted_shared_soln(s_s,:)
            end     
            '''
        elif opt_island['rep_pol'][num_isl] == 'random_n':           # Randomly Choose N
            '''
            available = 1:length(sorted_f)
            for s_s = 1:length(sorted_f)
                if s_s > OPT_island.rep_opt(num_isl)
                    break
                end
                if s_s > Npop
                    break
                end
                choice = available(random_num(1,length(available),'int'))
                popn(s_s,:) = sorted_shared_soln(choice,:)
                available(choice) = []
            end 
            '''
        elif opt_island['rep_pol'][num_isl] == 'best_n':             # Best N
           
            for s_s in range(opt_island['rep_opt'][num_isl]):
                if s_s > Npop:
                    break
                popn[s_s,:] = sorted_shared_soln[s_s,:]
           
        elif opt_island['rep_pol'][num_isl] == 'threshold_cost':     # Threshold Cost
            '''
            s_s = 1
            while sorted_f(s_s) <= OPT_island.rep_opt(num_isl) && s_s <= Npop
                popn(s_s,:) = sorted_shared_soln(s_s,:)
                s_s = s_s + 1
                if s_s > length(sorted_f)
                    break
                end
            end 
            '''
        elif opt_island['rep_pol'][num_isl] == 'threshold_percent':  # Threshold Percent (0-1)
            '''
            Ntake = round(OPT_island.rel_opt(num_isl)*length(sorted_f))
            for s_s = 1:Ntake
                if s_s > Npop
                    break
                end
                if s_s > length(sorted_f)
                    break
                end
                popn(s_s,:) = sorted_shared_soln(s_s,:) 
            end
            '''    
        else:
           
            raise Exception('Incorrect replacement method selected.\n')
            return
               
    return popn