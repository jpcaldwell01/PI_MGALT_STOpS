# FORM: var = MGALT_bestLegs(var,popn,f_trans)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function keeps track of the Best-So-Far Legs found
# |      throughout optimization
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -var               (10)       [dict]        [unitless]
# |         A dict containing the variable limits
# |     -f_trans           (Npop,transfers) [array]        [unitless]
# |         Evaluated cost of individual transfer legs for each popn member
# |     -popn              (Npop,len(member))       [int]           [unitless]
# |         The current migration's population
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -var          (10)       [struct]        [unitless]
# |         A dict containing the variable limits
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import numpy as np

def MGALT_bestLegs(var,popn,f_trans):
    
    for it1 in range(var['transfers']):
        for it2 in range(np.shape(f_trans)[0]):
            
            # if best leg arrays haven't been filled, append popn and costs. If these arrays were initialized with zeros, there is a possibility that bogus values could be substituted into optimization
            if np.shape(var['bestLegCosts']['leg%i' %(it1)])[0] < var['bestLegsKept']:
                var['bestLegCosts']['leg%i' %(it1)] = np.append(var['bestLegCosts']['leg%i' %(it1)],f_trans[it2,:][None,:],axis=0)
                var['bestLegMembers']['leg%i' %(it1)] = np.append(var['bestLegMembers']['leg%i' %(it1)],popn[it2,:][None,:],axis=0)
                continue
            # if popn[it2,I[it1]:I[it1]] in var['bestLegMembers']['leg%i' %(it1)] or popn[it2,I[0]:] in var['bestLegMembers']['leg%i' %(it1)]:
            max_f = np.max(var['bestLegCosts']['leg%i' %(it1)][:,it1])
            if f_trans[it2,it1] < max_f:
                max_f_I = np.argmax(var['bestLegCosts']['leg%i' %(it1)][:,it1])
                var['bestLegCosts']['leg%i' %(it1)][max_f_I,:]= f_trans[it2,:]
                var['bestLegMembers']['leg%i' %(it1)][max_f_I,:]= popn[it2,:]
                
    return var