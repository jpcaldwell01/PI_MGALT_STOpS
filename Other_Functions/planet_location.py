# %%

from numpy import subtract
from Other_Functions.constants import constants
from jplephem.spk import SPK

def planetLocation( planet,JD,opts ):
# # # LOOK UP PLANETARY STATE VECTORS # # #
# =============================================================================
# INPUTS:
#       planet:     planet ID [int or string]
#       JD:         time planet location is desired [JD]
# OUPUTS:
#   ~ return pos, vel ~
#       pos:        planet position wrt sun [km]
#       vel:        planet velocity wrt sun [km/s]
# =============================================================================
       
    # pre-allocate JPL data
    PlEph = SPK.open('de432s.bsp') # initiate JPL look-up

    # sun position and velocity wrt SS barycenter at JD
    posSun, velSun = PlEph[0,10].compute_and_differentiate(JD)

    # planet position and velocity wrt barycenter at JD
    pos, vel = PlEph[0,planet].compute_and_differentiate(JD) # [km], [km/day]

    # obtain planet location wrt sun, not SSB
    pos = subtract(pos,posSun)      # [km wrt. sun]

    # obtain planetary velociy wrt sun, not SSB
    vel = subtract(vel,velSun)      # [km/day]
    vel = vel/(constants('s2day')[0])  # [km/s]
    return pos, vel
