#  FORM: opts['PSO'] = parametersPSO(selection)
# |-----------------------------------------------------------------------
# | NOTES:
# |     -This function call outputs the necessary struct object to fully
# |     define a Particle Swarm Optimization Island. A switch/case allows 
# |     the user to add in commonly used cases to have it as an easy 
# |     return, and more cases can easily be added. The function also 
# |     allows for custom inputs, which is then exported to be in the 
# |     correct form for the Island struct object.
# |
# |     -A set of parameters is located below in the section titled 'MISC'
# |     which describes the object struct and the different parts of the 
# |     object and what they do. For more details on the methods and 
# |     inputs here the user should refer to Appendix A of the 
# |     accompanying thesis.
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -selection          (1,n)       [string]        [unitless]
# |         A used defined string which is used in a switch/case format
# |         to return predefined options
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -options            (1,n)       [assorted]      [unitless]
# |         The return struct with all the associated parameters
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |     The parameters are defined as follows:
# |     -Npop:              (1,1)       [int]           [unitless]
# |         Number of bees (solutions) in each population.
# |     -vmax:              (1,1)       [float]         [unitless]
# |         The max velocity for each variable. This number gets 
# |         multiplied by the max/min range for each variable to create 
# |         a max speed specific to each variables. The velocity is on a 
# |         scale from 0-1 but it is recommended to be at least 0.5 so 
# |         the particles start out moving at least halfway across the 
# |         variable space. 
# |     -tspan:             (1,1)       [int]           [unitless]
# |         How many time iterations are evaluated. This is simliar to 
# |         the number of gerations for the the Genetic Algorithm and 
# |         Differential Evolution.
# |     -K:                 (1,1)       [int]           [unitless]
# |         Number of informants. This is the number of bees that share 
# |         their best solution with all the other bees. It will choose 
# |         the best K solutions from the whole population.
# |     -cl:                (1,1)       [float]         [unitless]
# |         The bees' condfidence in their own velocity. A higher number 
# |         will encourage the bees to explore the search space. This 
# |         variables is on a 0-1 scale. Note that a bee's velocity will 
# |         only decrease over time.
# |     -cmax:              (1,1)       [float]         [unitless]
# |         The bees' confidence in other best solutions from the 
# |         informants. If this is too high it make cause premature 
# |         convergence. This variables is also on a 0-1 scale.
# |
# |-----------------------------------------------------------------------


# %%PSO Parameters Selection
def parametersPSO(selection):
    options = {}
    if selection == '50_50':     #50 popn and 50 iter

        options['Npop'] = 50    # Number of members in each population
        options['tspan'] = 50   # How many iterations to run for
        options['vmax'] = 0.7   # The max velocity of each particle
        options['K'] = 5        # The number of informants
        options['cl'] = 0.9     # The confidence in their own velocity
        options['cmax'] = 0.7   # The confidence in other solutions

    elif selection == '50_75':   #50 popn and 75 iter

        options['Npop'] = 50    # Number of members in each population
        options['tspan'] = 75   # How many iterations to run for
        options['vmax'] = 0.7   # The max velocity of each particle
        options['K'] = 4        # The number of informants
        options['cl'] = 0.9     # The confidence in their own velocity
        options['cmax'] = 0.7   # The confidence in other solutions
        
    elif selection == '250_100': #250 popn and 100 iter

        options['Npop'] = 250    # Number of members in each population
        options['tspan'] = 100   # How many iterations to run for
        options['vmax'] = 0.7   # The max velocity of each particle
        options['K'] = 5        # The number of informants
        options['cl'] = 0.9     # The confidence in their own velocity
        options['cmax'] = 0.7   # The confidence in other solutions
    
    elif selection == '350_100': #350 popn and 100 iter

        options['Npop'] = 350    # Number of members in each population
        options['tspan'] = 100   # How many iterations to run for
        options['vmax'] = 0.7   # The max velocity of each particle
        options['K'] = 4        # The number of informants
        options['cl'] = 0.9     # The confidence in their own velocity
        options['cmax'] = 0.7   # The confidence in other solutions
    
    elif selection == '350_150': #350 popn and 150 iter

        options['Npop'] = 350    # Number of members in each population
        options['tspan'] = 150   # How many iterations to run for
        options['vmax'] = 0.7   # The max velocity of each particle
        options['K'] = 5        # The number of informants
        options['cl'] = 0.9     # The confidence in their own velocity
        options['cmax'] = 0.7   # The confidence in other solutions
    else:
        raise Exception("Ivalid PSO Parameter")
                
    return options
