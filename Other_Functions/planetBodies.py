#  |-----------------------------------------------------------------------
#  |
#  | NOTES:
#  |     -Return R, V, and JD arrays for bodies included in the problem for 
#  |     the valid time span
#  |
#  |-----------------------------------------------------------------------
#  |
#  | INPUTS:
#  |     -BOD                (1,1)       [struct]        [unitless]
#  |         A struct containing information pertaining to the planetary
#  |         bodies. Contains list of bodies, launch windows and ToF, and 
#  |         planetary R/V/JD vectors. This struct has dynamic fields and 
#  |         will adapt to contain only the necesary information
#  |
#  |-----------------------------------------------------------------------
#  |
#  | OUTPUTS:
#  |     -bodies_R     	(3*Nbodies,n)	[float]         [AU]
#  |         Heliocentric Radius components for all planets
#  |     -bodies_V       (3*Nbodies,n)	[float]         [AU/TU]
#  |         Heliocentric Velocity components for all planets
#  |     -bodies_JD   	(Nbodies,n)     [float]         [JD]
#  |         Julian day for the departure body and all transfer bodies
#  |
#  |-----------------------------------------------------------------------
#  |
#  | MISC:
#  |
#  |-----------------------------------------------------------------------

import datetime
from julian import to_jd, from_jd
from Other_Functions.planet_location import planetLocation
import numpy as np
import calendar
import scipy.io

# %% Initial
def planetBodies(bod,opts):
    bodies = bod['bodies']
    bodiesInt = [0 for its in range(len(bodies))]
    for it in range(len(bodies)):
        if bodies[it] == 'Sun':
            bodiesInt[it] = 0
        elif bodies[it] == 'Mercury':
            bodiesInt[it] = 1
        elif bodies[it] == 'Venus':
            bodiesInt[it] = 2       
        elif bodies[it] == 'Earth':
            bodiesInt[it] = 3               
        elif bodies[it] == 'Mars':
            bodiesInt[it] = 4        
        elif bodies[it] == 'Jupiter':
            bodiesInt[it] = 5        
        elif bodies[it] == 'Saturn':
            bodiesInt[it] = 6        
        elif bodies[it] == 'Uranus':
            bodiesInt[it] = 7        
        elif bodies[it] == 'Neptune':
            bodiesInt[it] = 8        
        elif bodies[it] == 'Pluto':
            bodiesInt[it] = 9        
        
    tspan = bod['tspan']
    
    
    # Calculate Start/End year based off tspan
    date1 = from_jd(tspan[0])
    date2 = from_jd(tspan[1])
    year1 = date1.year
    year2 = date2.year
    year_span = year2-year1
    
    # Add the years between the start and ending year
    if year_span > 1:    # If larger than 1 year
        
        year = [0]*(year2-year1+1)
        
        for i1 in range(len(year)):
            year[i1] = year1+i1
        
    else:
        year = [year1]
    
    # Preallocate
    daysNum = [0 for i in range(len(year))]
    Data = [[{'JD','Date','R','V'}for i in range(len(bodies))] for i in range(len(year))]
    
    for i in range(len(year)):
        daysNum[i] = 366 if calendar.isleap(year[i]) else 365

    total_indices = 0
    array_indices = [0,0]
    
    # Load Appropriate Data
    for i2 in range(len(year)):
        JD = [0 for it in range(daysNum[i2]*4)]
        date = [0 for it in range(len(JD))]
        it = 0
        for it in range(daysNum[i2]*4):
            JD[it] = to_jd(datetime.datetime(year[i2],1,1)) + it/4
            date[it] = from_jd(JD[it])

        for i3 in range(len(bodies)):
            
            if opts['ephemeris'] == 'jplEphem':
                
                r,v = planetLocation(bodiesInt[i3],JD,opts)
                Data[i2][i3] = {'JD' : JD, 'Date' : date,'R':r,'V':v}
                
            elif opts['ephemeris'] == 'JPL_Horizons':
                
                string = 'PlanetData/'+bodies[i3]+'_'+str(year[i2])+'.mat'
                tempStuff = scipy.io.loadmat(string)
                Data[i2][i3] = {'JD' : tempStuff['JD'].T, 'Date' : tempStuff['DATE'],'R':tempStuff['R'],'V':tempStuff['V']}
            
            # Add something for when data is not available (ode45)
            
    # Make a array of the matrix indicies for output data
    total_indices = sum(daysNum)*4
    
    # Pre allocate

    bodies_JD = np.array([])
    bodies_R = np.empty((len(bodies)*3,0))
    bodies_V = bodies_R
    
    # Put the data from the extracted files into the return values

    for i5 in range(len(year)):
        JDtemp = np.array(Data[i5][0]['JD'])
        bodies_JD = np.append(bodies_JD,JDtemp)
        
        Rtemp = np.empty([len(bodies)*3,len(JDtemp)])
        Vtemp = np.empty([len(bodies)*3,len(JDtemp)])
        for i6 in range(len(bodies)):
            Rtemp[(i6+1)*3-3:(i6+1)*3,:] = np.array(Data[i5][i6]['R'])
            Vtemp[(i6+1)*3-3:(i6+1)*3,:] = np.array(Data[i5][i6]['V'])
    
        bodies_R = np.append(bodies_R,Rtemp,axis=1)
        bodies_V = np.append(bodies_V,Vtemp,axis=1)        
    bodies_JD = bodies_JD.tolist()
    bodies_R = bodies_R.tolist()
    bodies_V = bodies_V.tolist()
    
    return bodies_R, bodies_V, bodies_JD