# FORM: [const] = getconstants(BOD)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -Get the constants for this problem
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -BOD                (1,1)       [struct]        [unitless]
# |         A struct containing information pertaining to the planetary
# |         bodies. Contains list of bodies, launch windows and ToF, and 
# |         planetary R/V/JD vectors. This struct has dynamic fields and 
# |         will adapt to contain only the necesary information
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -const              (1,1)       [struct]        [unitless]
# |         A struct containing constants used in the calcs. Contains
# |         values for AU, TU, Sun (rad/mu/rp) and (rad/mu/rp/SOI/per) 
# |         for any bodies used in the optimization scheme. This is a 
# |         dynamic struct and will adapt to contain only the necesary 
# |         information
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

from Other_Functions.constants import constants
def getConstants(bod):
    ## Get constants
    
    bodies = bod['bodies']
    const = {}
    
    # Get the canonical units
    const['AU'] = constants('AU')
    const['TU'] = constants('TU')
    
    # Get sun units
    const['Sun_rad'] = constants('radius Sun')
    const['Sun_mu'] = constants('mu Sun')
    const['Sun_rp'] = constants('min rp Sun')
    
    # Get all related info to bodies
    for i1 in range(len(bodies)):
            
        if bodies[i1] == 'Mercury':
            
            const['Mercury_rad'] = constants('radius Mercury')
            const['Mercury_mu'] = constants('mu Mercury')
            const['Mercury_rp'] = constants('min rp Mercury')
            const['Mercury_SOI'] = constants('SOI Mercury')
            const['Mercury_per'] = constants('period Mercury')
            
        elif bodies[i1] == 'Venus':
            
            const['Venus_rad'] = constants('radius Venus')
            const['Venus_mu'] = constants('mu Venus')
            const['Venus_rp'] = constants('min rp Venus')
            const['Venus_SOI'] = constants('SOI Venus')
            const['Venus_per'] = constants('period Venus')
            
        elif bodies[i1] == 'Earth':
            
            const['Earth_rad'] = constants('radius Earth')
            const['Earth_mu'] = constants('mu Earth')
            const['Earth_rp'] = constants('min rp Earth')
            const['Earth_SOI'] = constants('SOI Earth')
            const['Earth_per'] = constants('period Earth')
            
        elif bodies[i1] == 'Mars':
            
            const['Mars_rad'] = constants('radius Mars')
            const['Mars_mu'] = constants('mu Mars')
            const['Mars_rp'] = constants('min rp Mars')
            const['Mars_SOI'] = constants('SOI Mars')
            const['Mars_per'] = constants('period Mars')
            
        elif bodies[i1] == 'Jupiter':
            
            const['Jupiter_rad'] = constants('radius Jupiter')
            const['Jupiter_mu'] = constants('mu Jupiter')
            const['Jupiter_rp'] = constants('min rp Jupiter')
            const['Jupiter_SOI'] = constants('SOI Jupiter')
            const['Jupiter_per'] = constants('period Jupiter')
            
        elif bodies[i1] == 'Saturn':
            
            const['Saturn_rad'] = constants('radius Saturn')
            const['Saturn_mu'] = constants('mu Saturn')
            const['Saturn_rp'] = constants('min rp Saturn')
            const['Saturn_SOI'] = constants('SOI Saturn')
            const['Saturn_per'] = constants('period Saturn')
            
        elif bodies[i1] == 'Uranus':
            
            const['Uranus_rad'] = constants('radius Uranus')
            const['Uranus_mu'] = constants('mu Uranus')
            const['Uranus_rp'] = constants('min rp Uranus')
            const['Uranus_SOI'] = constants('SOI Uranus')
            const['Uranus_per'] = constants('period Uranus')
            
        elif bodies[i1] == 'Neptune':
            
            const['Neptune_rad'] = constants('radius Neptune')
            const['Neptune_mu'] = constants('mu Neptune')
            const['Neptune_rp'] = constants('min rp Neptune')
            const['Neptune_SOI'] = constants('SOI Neptune')
            const['Neptune_per'] = constants('period Neptune')
            
        elif bodies[i1] == 'Pluto':
            
            const['Pluto_rad'] = constants('radius Pluto')
            const['Pluto_mu'] = constants('mu Pluto')
            const['Pluto_rp'] = constants('min rp Pluto')
            const['Pluto_SOI'] = constants('SOI Pluto')
            const['Pluto_per'] = constants('period Pluto')
            
        else:
            raise Exception('Incorrect body selection.')
            return
    
    return const


