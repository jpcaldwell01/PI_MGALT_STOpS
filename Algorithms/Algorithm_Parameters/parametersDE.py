#  FORM: opts['DE'] = parametersDE(selection)
# |-----------------------------------------------------------------------
# | NOTES:
# |     -This function call outputs the necessary struct object to fully
# |     define a Differential Evolution Island. A switch/case allows the 
# |     user to add in commonly used cases to have it as an easy return, 
# |     and more cases can easily be added. The function also allows for 
# |     custom inputs, which is then exported to be in the correct form 
# |     for the Island struct object.
# |
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
# |         Number of members (solutions) in each population.
# |     -Ngen:              (1,1)       [int]           [unitless]
# |         Number of generations that will be evaluated.
# |     -pc:                (1,1)       [float]         [percent]
# |         Percent probability for a crossover to occur on a 0-1 scale.
# |     -sel_method:        (1,1)       [string]        [unitless]
# |         Selection method for the base vector. Options are: 'random', 
# |         'best_so_far', 'random_best_blend'. All options are case 
# |         sensitive.
# |     -F_method:          (1,1)       [string]        [unitless]
# |         Method for applying the scale factor. Options are: 'constant', 
# |         'jitter', and 'dither'. All options are case sensitive.
# |     -F:                 (1,1/2)     [float]         [unitless]
# |         Scaling factor when applying difference vector to base vector. 
# |         If  F  method is 'constant' then this is a (1,1) from 0-1. 
# |         For F methods of 'jitter' and 'dither' the F value is a range 
# |         so the input is (1,2) with the first input being the low end 
# |         of the range and the second input being the high end of the 
# |         range. These values are also 0-1. It is recommended in both 
# |         cases that F is at least 0.4.
# |     -surv_method:       (1,1)       [string]        [unitless]
# |         Method for choosing survivors. Options are: 
# |         'natural_selection' and 'tournament'. All options are case 
# |         sensitive.
# |     -T:                 (1,1)       [int]           [unitless]
# |         The number of competitors in each tournament if the survivor 
# |         method is 'tournament'. Must be 2 or greater.
# |     -weight:            (1,1)       [string]        [unitless]
# |         How to calculate selection probablilities. Options are: 'cost' 
# |         and 'rank'. All options are case sensitive.
# |
# |-----------------------------------------------------------------------

# %% DE Parameters Selection
def parametersDE(selection):
    options = {}

    if selection == '75_30':      #75 popn and 30 gen
        options['Npop'] = 75                     	   # Number of members in each population
        options['Ngen'] = 30                           # Number of generations to evaluate
        options['pc'] = 0.8                            # Percent probability for a crossover to occur
        options['sel_method'] = 'random_best_blend'    # The selection method for the base vector
        options['F_method'] = 'jitter'                 # Method for applying the scale factor
        options['F'] = (0.5,0.9)                       # The scaling factor between vector and base vector
        options['surv_method'] = 'tournament'          # The survivor method         
        options['T'] = 3                               # Number of competitors in each tourneyment
        options['weight'] = 'cost'                     # How to calculate selection probabilities
        
    elif selection == '100_30':     #100 popn and 30 gen

        options['Npop'] = 100                     	   # Number of members in each population
        options['Ngen'] = 30                           # Number of generations to evaluate
        options['pc'] = 0.8                            # Percent probability for a crossover to occur
        options['sel_method'] = 'random_best_blend'    # The selection method for the base vector
        options['F_method'] = 'jitter'                 # Method for applying the scale factor
        options['F'] = (0.5,0.9)                       # The scaling factor between vector and base vector
        options['surv_method'] = 'tournament'          # The survivor method         
        options['T'] = 2                               # Number of competitors in each tourneyment
        options['weight'] = 'cost'                     # How to calculate selection probabilities
        
    elif selection == '300_115':     #350 popn and 75 gen
        
        options['Npop'] = 300                     	   # Number of members in each population
        options['Ngen'] = 115                          # Number of generations to evaluate
        options['pc'] = 0.8                            # Percent probability for a crossover to occur
        options['sel_method'] = 'random_best_blend'    # The selection method for the base vector
        options['F_method'] = 'jitter'                 # Method for applying the scale factor
        options['F'] = (0.5,0.9)                       # The scaling factor between vector and base vector
        options['surv_method'] = 'tournament'          # The survivor method         
        options['T'] = 3                               # Number of competitors in each tourneyment
        options['weight'] = 'cost'                     # How to calculate selection probabilities
        
    elif selection == '350_75':     #350 popn and 75 gen
        
        options['Npop'] = 350                     	   # Number of members in each population
        options['Ngen'] = 75                           # Number of generations to evaluate
        options['pc'] = 0.8                            # Percent probability for a crossover to occur
        options['sel_method'] = 'random_best_blend'    # The selection method for the base vector
        options['F_method'] = 'jitter'                 # Method for applying the scale factor
        options['F'] = (0.5,0.9)                       # The scaling factor between vector and base vector
        options['surv_method'] = 'tournament'          # The survivor method         
        options['T'] = 2                               # Number of competitors in each tourneyment
        options['weight'] = 'cost'                     # How to calculate selection probabilities
    
    elif selection == '450_180':     #450 popn and 180 gen

       options['Npop'] = 450                     	   # Number of members in each population
       options['Ngen'] = 180                           # Number of generations to evaluate
       options['pc'] = 0.8                            # Percent probability for a crossover to occur
       options['sel_method'] = 'random_best_blend'    # The selection method for the base vector
       options['F_method'] = 'jitter'                 # Method for applying the scale factor
       options['F'] = (0.5,0.9)                       # The scaling factor between vector and base vector
       options['surv_method'] = 'tournament'          # The survivor method         
       options['T'] = 3                              # Number of competitors in each tourneyment
       options['weight'] = 'cost'                     # How to calculate selection probabilities
    else:
        raise Exception("Invalid DE Parameter")
        
    return options

