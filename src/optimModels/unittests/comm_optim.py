import os
import optimModels
from framed import load_cbmodel
from optimModels.comm_optim.CModel import CModel
from optimModels.comm_optim.ea_setup import change_config

if __name__ == "__main__":
    # Step 1
    # Load any ammount of models in framed objects
    # All of the models should be loaded in the same notation style
    # The 10 models loaded below were arbitrarily chosen
    optimmodels_path = os.path.dirname(optimModels.__file__)
    models_path = os.path.abspath(os.path.join(optimmodels_path, os.pardir, os.pardir, "examples", "models"))

    model1 = "Yokenella_regensburgei_ATCC_43003.xml"
    model2 = "Acinetobacter_junii_SH205.xml"
    model3 = "Clostridiales_sp_1_7_47FAA.xml"
    model4 = "Achromobacter_xylosoxidans_A8.xml"
    model5 = "Achromobacter_xylosoxidans_NBRC_15126.xml"
    model6 = "Acidaminococcus_intestini_RyC_MR95.xml"
    model7 = "Acidaminococcus_sp_D21.xml"
    model8 = "Acinetobacter_calcoaceticus_PHEA_2.xml"
    model9 = "Acinetobacter_lwoffii_WJ10621.xml"
    model10 = "Actinobacillus_pleuropneumoniae_L20.xml"
    model_list = [model1, model2, model3, model4, model5, model6, model7, model8, model9, model10]

    list_models = [load_cbmodel(filename = str(os.path.join(models_path, model)),
                                flavor = "cobra:other") for model in model_list]

    # Step 2
    # Create the Community Model
    comm_model = CModel(
        community_id = "model_id",
        models = list_models,
        empty_flag = False  # This creates a complete medium, to crete an empty medium input True
        )                                       # That way you can open any medium through constraints

    # Step 3
    # Run the EA Optimization
    # 3.1 Configure the optimization parameters
    change_config(  # Leave parameters as None to use defaults
        pop_size = None,  # changes the size of the populations of candidates
        max_gen = None,  # changes the max number of generations
        mut = None,  # changes the mutation rate when creating new solutions, 0 < x < 1
        cross = None  # changes the crossover rate when creating new solutions, 0 < x < 1
        )

    # 3.2 Configure the type of optimization
    ea_optim_options = [  # this is a list with two tuples
        (  # this first tuple configures the inital population
            "random",  # "random" is random, "headstart" forces  some organisms in all solutions
            [1, 0, 1, 0, 0, 0, 0, 0, 0, 0]  # only with headstart, needs sample, see comment after the list*
            ),
        (  # this tuple conifgures the creation of the follwing generations
            "tourn",  # "tourn" is default, "changeworst" replaces the worst, "keep" pairs with "headstart"
            [1, 0, 1, 0, 0, 0, 0, 0, 0, 0]  # changeworst needs integer, keep needs sample, see comment after*
            )
        ]

    # *the second position of each tuple is used for some of the options in the first position of each tuple
    # when selecting headstart or keep, you need to provide a sample of the organisms you want to keep
    # the format of this sample, however, changes with the quantity option in the ea function below
    # if the quantity is set at 0, then the sample is an integer list as set above. this is list should be the same
    # size of the number of models in the CModel object. in that list you should place 1 in the positions of the
    # models you want to keep, and 0 in the positions that are to be changeable.
    # if, however, you set a integer value to the quantity in the ea function, the sample is in integer format
    # intead of binary, where you select, directly the positions you want to keep
    # meaning, [1, 0, 1, 0, 0, 0, 0, 0, 0, 0] = [0, 2]
    # when using changeworst, to create new candidates, the algorithm takes exhisting ones, checks,
    # the performance of each individual organism and swaps out an x number of the worst performing
    # organisms for the objective provided. This x is the integer you should place in the second position
    # of the second tuple in that case in the next lines some examples will be commented

    options_1 = [("random", None), ("tourn", None)]  # default settings, completely random
    options_2a = [("headstart", [1, 0, 1, 0, 0, 0, 0, 0, 0, 0]), ("keep", [1, 0, 1, 0, 0, 0, 0, 0, 0, 0])]  # keeps
    # organisms  in all solutions, when quantity = 0
    options_2b = [("headstart", [0, 2]), ("keep", [0, 2])]  # equal to 2a but for quantity > 0
    options_3 = [("random", None), ("changeworst", 2)]  # changes the 2 worst performing organisms of
    # each candidate every generation

    # 3.3 Run the Optimization
    ea_optim_result = comm_model.ea(
        options = options_1,
        quantity = 3,
        obj_list = [],  # list of objective reactions str ids
        cons = {},  # dictionary of constraints, format {"reac_id" : (lower bound, upper bound), ...}
        fitness = None,  # reactino id str to use as fitness instead of the objective
        goal = None  # goal to stop algorithm before the max gens if desired value is reached
        )

    print("Best Solution: ", ea_optim_result)
