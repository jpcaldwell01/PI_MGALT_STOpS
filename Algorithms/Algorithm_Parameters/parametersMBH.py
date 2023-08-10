#  FORM: opts['MBH'] = parametersMBH(selection)
# |-----------------------------------------------------------------------
# | NOTES:
# |     -This function call outputs the necessary struct object to fully
# |     define a Monotonic Basin Hopping Island. A switch/case allows the user 
# |     to add in commonly used cases to have it as an easy return, and 
# |     more cases can easily be added. The function also allows for 
# |     custom inputs, which is then exported to be in the correct form 
# |     for the Island struct object.
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
# |     -N1_Outer          (1,1)       [int]           [unitless]
# |         Number of members in each initial population.
# |     -N1_Inner          (1,1)       [int]           [unitless]
# |         Number of iterations to run in initial local search.
# |     -N2_Outer          (1,1)       [int]           [unitless]
# |         Number of members to generate from each parent member for basin exploration.
# |     -N2_Inner           (1,1)       [int]           [unitless]
# |         Number of members to generate from each parent member for basin exploration.
# |     -maxclst            (1,1)       [int]           [unitless]
# |         Maximum number of clusters to explore if ‘N1_Outer’ deemed feasible
# |     -per_feas           (1,Nmig+1)     [Array of Float]       [%]
# |         Percentage values to apply to the parameter used to check feasibility.
# |     -per_rand           (1,Nmig+1)     [Array of Float]       [%]
# |         Percentage values to apply to the perturbation limits generation of new members
# |     
# |
# |-----------------------------------------------------------------------


def parametersMBH(selection):
   
    options = {}    # MBH Parameters

    if selection == '15000_1500':
        options['N1_Outer'] = 15000         # Number of members in each initial population
        options['N1_Inner'] = 1500          # Number of iterations to run in initial local search.
        options['N2_Outer'] = 200           # Number of members to generate from each parent member for basin exploration.
        options['N2_Inner'] = 125           # Number of members to generate from each parent member for basin exploration.
        options['maxclst'] = 5              # Maximum number of clusters to explore if ‘N1 outer’ deemed feasible
        options['per_feas'] = [75.0,10.0,2.20,1.30,1.00] # Percentage values to apply to the parameter used to check feasibility.
        options['per_rand'] = [0.13,0.09,0.06,0.03,0.01] # Percentage values to apply to the perturbation limits generation of new members
        options['feas_check'] = any        #
        options['feas_tol'] = 1e-5         #
    
    elif selection == '1450_75':
        options['N1_Outer'] = 1450          # Number of members in each initial population
        options['N1_Inner'] = 75            # Number of iterations to run in initial local search.
        options['N2_Outer'] = 175           # Number of members to generate from each parent member for basin exploration.
        options['N2_Inner'] = 100           # Number of members to generate from each parent member for basin exploration.
        options['maxclst'] = 5              # Maximum number of clusters to explore if ‘N1 outer’ deemed feasible
        options['per_feas'] = [75.0,10.0,2.20] # Percentage values to apply to the parameter used to check feasibility.
        options['per_rand'] = [0.13,0.09,0.06] # Percentage values to apply to the perturbation limits generation of new members
        options['feas_check'] = any        #
        options['feas_tol'] = 1e-5         #
        
    elif selection == '2000_350':
        options['N1_Outer'] = 2000          # Number of members in each initial population
        options['N1_Inner'] = 350            # Number of iterations to run in initial local search.
        options['N2_Outer'] = 100           # Number of members to generate from each parent member for basin exploration.
        options['N2_Inner'] = 25           # Number of members to generate from each parent member for basin exploration.
        options['maxclst'] = 5              # Maximum number of clusters to explore if ‘N1 outer’ deemed feasible
        options['per_feas'] = [5.00,2.50,1.75,1.25,1.00] # Percentage values to apply to the parameter used to check feasibility.
        options['per_rand'] = [.1,.07,.05,.02,.01] # Percentage values to apply to the perturbation limits generation of new members
        options['feas_check'] = any        #
        options['feas_tol'] = 1e-5         #
    
    elif selection == '750_300':
        options['N1_Outer'] = 750           # Number of members in each initial population
        options['N1_Inner'] = 300           # Number of iterations to run in initial local search.
        options['N2_Outer'] = 175           # Number of members to generate from each parent member for basin exploration.
        options['N2_Inner'] = 100           # Number of members to generate from each parent member for basin exploration.
        options['maxclst'] = 5              # Maximum number of clusters to explore if ‘N1 outer’ deemed feasible
        options['per_feas'] = [75.0,10.0,2.20] # Percentage values to apply to the parameter used to check feasibility.
        options['per_rand'] = [0.13,0.09,0.06] # Percentage values to apply to the perturbation limits generation of new members
        options['feas_check'] = any        #
        options['feas_tol'] = 1e-5         #
        
    elif selection == '10000_800':
        options['N1_Outer'] = 10000         # Number of members in each initial population
        options['N1_Inner'] = 800           # Number of iterations to run in initial local search.
        options['N2_Outer'] = 175           # Number of members to generate from each parent member for basin exploration.
        options['N2_Inner'] = 100           # Number of members to generate from each parent member for basin exploration.
        options['maxclst'] = 5              # Maximum number of clusters to explore if ‘N1 outer’ deemed feasible
        options['per_feas'] = [75.0,10.0,2.20,1.30]#,1.00] # Percentage values to apply to the parameter used to check feasibility.
        options['per_rand'] = [0.13,0.09,0.06,0.03]#,0.01] # Percentage values to apply to the perturbation limits generation of new members
        options['feas_check'] = any        #
        options['feas_tol'] = 1e-5         #
        
    else:
        print(selection,' is not a valid preset for MBH Parameters. Please try again.')
        raise Exception()
    
    return options