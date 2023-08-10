MGALT-STOpS 1.3
===============

MGALT-STOpS is a software suite which allows for optimization of low-thrust spacecraft trajectories between multiple planets.
The Spacecraft Trajectory Optimization Suite (STOpS) uses Multiple Gravity-Assist Low-Thrust (MGALT) trajectories paired with the island model paradigm to accomplish this goal. 
The island model utilizes four different global search algorithms: a Genetic Algorithm, Differential Evolution, Particle Swarm Optimization, and Monotonic Basin Hopping.
Each island runs either a direct or indirect optimization method to solve for the most optimal trajectory between multiple planets.


## Examples

Included in this download is a demo program, titled "PI_MGALT_STOpS_DEMO", which provides a step by step walkthrough of MGALT-STOpS and explains everything along the way.
It is recommended to become familiar with the demo script before venturing into customization with the main script.


## Framework

All of the files were originally built and tested with MATLAB 2019b [MATLAB](https://www.mathworks.com/products/matlab.html).
New to Version 1.3, MGALT-STOpS has been mostly converted to Python 3, specifically for use with the Spider IDE through Anaconda (https://www.anaconda.com/download)

There are two modules that need to be downloaded: julian (https://anaconda.org/conda-forge/julian), and 

## This Version

This version ported over most of the MATLAB suite into Python 3. There reamins some code left to port, particularly the LT_DIR_FSM_2D and LT_IN_FSM_2D solvers and related functions.
This version also included various pruning features, allowing users to reduce their search space reliably without risk of pruning away good solutions.

## More Information

This code was inherited/developed/updated by Caldwell to support his thesis, "Performance Improvements for the Multiple Gravity Assist Low Thrust Spacecraft Trajectory Optimization Suite (PI-MGALT-STOpS)"; a requirement for a MS in Aerospace Engineering at California Polytechnic State University, San Luis Obispo.

The development of MGALT-STOpS was a continuation of prior work on STOpS by Shane P. Sheehan and Timothy J. Fitzgerald.
The supporting thesis for PI-MGALT-STOpS, MGALT-STOpS, Low-Thrust STOpS, and STOpS can be found below:
* ["Performance Improvements for the Multiple Gravity Assist Low Thrust Spacecraft Trajectory Optimization Suite (PI-MGALT-STOpS)"](Publication Pending, pdf located in repository) by Caldwell (2023)
* ["Spacecraft Trajectory Optimization Suite (STOpS): Design and Optimization of Multiple Gravity-Assist Low-Thrust (MGALT) Trajectories Using Modern Optimization Techniques"](https://digitalcommons.calpoly.edu/theses/2247/) by MALLOY (2020)
* ["Spacecraft Trajectory Optimization Suite (STOPS): Optimization of Low-Thrust Interplanetary Spacecraft Trajectories Using Modern Optimization Techniques"](https://digitalcommons.calpoly.edu/theses/1901/) by SHEEHAN (2017)
* ["Spacecraft Trajectory Optimization Suite (STOpS): Optimization of Multiple Gravity Assist Spacecraft Trajectories Using Modern Optimization Techniques"](https://digitalcommons.calpoly.edu/theses/1503/) by FITZGERALD (2015)

The MGALT-STOpS Version 1.2 author can be contacted on Discord at "HuntrTrakr#8835" if there are any questions.
The MGALT-STOpS Version 1.3 author can be contacted on Discord at "Jypsee#2339" if there are any questions.

