import os
import optimModels
from collections import OrderedDict
from optimModels.utils.constantes import optimType
from optimModels.optimization.run import kinetic_strain_optim
from optimModels.model.kineticModel import load_kinetic_model
from optimModels.utils.configurations import KineticConfigurations
from optimModels.simulation.simul_problems import KineticSimulationProblem
from optimModels.optimization.evaluation_functions import build_evaluation_function


def build_chassaganole_map():
    """
    auxiliary function for the chassagnole model parameter map
    :return dict: paramter map in the format {reaction id: [parameter id, ...], ...}
    """

    optimmodels_path = os.path.dirname(optimModels.__file__)
    example_chassagnole_path = os.path.join(optimmodels_path, "examples", "models", "chassagnole2002.xml")

    model1 = load_kinetic_model(filename = example_chassagnole_path, kmap = {})
    # model1 = load_kinetic_model(filename = models_path + "chassagnole2002.xml", kmap = {})
    parameters = list(model1.get_parameters())
    params = parameters[10:]

    param_map = {}
    for param in params:
        reaction = param.split('_')[0]
        if reaction not in param_map.keys():
            param_map[reaction] = [param]

    return param_map


def kinetic_optimization(model, optim_type, eval_func, **kwargs):
    """
    default template function for kinetic model optimization
    :param load_kinetic_model() model: kinetic model object
    :param str optim_type: "KO" - Knockouts or "UO" - Under/Over expression
    :param build_evaluation_function() eval_func: evaluating function
    :param kwargs: all of the optional arguments
    :return kinetic_strain_optim(): result of the optimization
    """
    output_file = kwargs.get("output_file", "results.csv")
    time = kwargs.get("time", KineticConfigurations.STEADY_STATE_TIME)
    cand_size = kwargs.get("size", None)
    multi_thread = kwargs.get("isMultiProc", False)
    crit_params = kwargs.get("critical_parameters", [])

    kinetic_uo_levels = kwargs.get(
        "levels", [2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5])

    if optim_type == "KO":
        optim_type = optimType.REACTION_KO
    elif optim_type == "UO":
        optim_type = optimType.REACTION_UO

    kinetic_simul_problem = KineticSimulationProblem(model, tSteps=[0, time])
    res = kinetic_strain_optim(
        simulProblem = kinetic_simul_problem,
        objFunc = eval_func,
        levels = kinetic_uo_levels,
        type = optim_type,
        isMultiProc = multi_thread,
        candidateSize = cand_size,
        criticalParameters = crit_params,
        resultFile = output_file)
    return res


if __name__ == "__main__":
    # First Step:
    # Loading the model
    optimmodels_path = os.path.dirname(optimModels.__file__)
    example_chassagnole_path = os.path.join(optimmodels_path, "examples", "models", "chassagnole2002.xml")
    example_chassagnole_kmap = build_chassaganole_map()
    ex_kinetic_model = load_kinetic_model(
        filename = example_chassagnole_path,
        kmap = OrderedDict(example_chassagnole_kmap)
        )

    # Second Step
    # Evaluation Function
    example_evaluation = build_evaluation_function("targetFlux", ["vPK"])

    # Third Step
    # Optimization
    results = kinetic_optimization(
        model = ex_kinetic_model,
        optim_type = "KO",  # KO or KO
        eval_func = example_evaluation,
        critical_parameters = [],  # list of paramters ids str to avoid in optmization
        isMultiProc = False,  # replace with True to use multi-processing
        size = 10,  # chose max number of alterations, None uses default
        output_file = "optimization_results.csv"
        )
    for result in results:
        result.print()
