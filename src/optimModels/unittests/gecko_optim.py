import os.path
import geckopy
from optimModels.utils.constantes import optimType
from optimModels.optimization.run import gecko_strain_optim
from optimModels.simulation.simul_problems import GeckoSimulationProblem
from optimModels.optimization.evaluation_functions import build_evaluation_function
from optimModels.unittests.gecko_simul import convert_mmol_to_g, loading_yeast_gecko, loading_any_gecko


def gecko_optimization(model, optim_type, eval_func, **kwargs):
    """
    This function is a default template for any gecko optimization

    :param GeckoModel model: GeckoModel object from geckopy
    :param str optim_type: "KO" - Knockouts or "UO" - Under/Over expression
    :param build_evaluation_function() eval_func: evaluating function
    :param kwargs: all of the optional arguments
    :return: None
    """
    gecko_constraints = kwargs.get("constraints", [])
    gecko_simul_problem = GeckoSimulationProblem(model, constraints = gecko_constraints)
    gecko_wt_fluxes = gecko_simul_problem.simulate()
    gecko_simul_problem.wt_concentrations = gecko_wt_fluxes.get_protein_concentrations()

    gecko_uo_levels = kwargs.get(
        "levels", [2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5])

    if optim_type == "KO":
        optim_type = optimType.PROTEIN_KO
    elif optim_type == "UO":
        optim_type = optimType.PROTEIN_UO

    critical_proteins = kwargs.get("critical_proteins", [])
    if critical_proteins == "auto":
        critical_proteins = gecko_simul_problem.find_essential_proteins()

    multi_thread = kwargs.get("isMultiProc", False)
    cand_size = kwargs.get("size", None)
    output_file = kwargs.get("output_file", "results.csv")

    res = gecko_strain_optim(
        simulProblem = gecko_simul_problem,
        evaluationFunc = eval_func,
        levels = gecko_uo_levels,
        type = optim_type,
        criticalProteins = critical_proteins,
        isMultiProc = multi_thread,
        candidateSize = cand_size,
        resultFile = output_file
        )
    return res


if __name__ == "__main__":
    # First Step:
    # Loading the model
    # Using the provided yeast gecko model:
    yeast_single_pool = loading_yeast_gecko()

    # Using the multi-pool variant:
    gecko_path = os.path.dirname(geckopy.__file__)
    protein_ggdw_path = os.path.join(gecko_path, "data_files", "sanchez-mmol_gdw.csv")
    protein_ggdw = convert_mmol_to_g(protein_ggdw_path)

    yeast_multi_pool = loading_yeast_gecko(prot_measure_ggdw = protein_ggdw)

    # Using an outside model:
    path_to_any_gecko = "fake_path/any_gecko.xml"
    any_gecko_biomass = "any_biomass"
    any_gecko_single_pool = loading_any_gecko(path = path_to_any_gecko,
                                              biomass = any_gecko_biomass)
    # turn into multi_pool using limit_proteins (also possible to use franctions)
    path_to_any_ggdw = "fake_path/any_ggdw.csv"
    any_protein_ggdw = convert_mmol_to_g(path_to_any_ggdw)
    any_gecko_multi_pool = any_gecko_single_pool.limit_proteins(ggdw = any_protein_ggdw)

    # Second step:
    # Evaluating Function         # different functions require different arguments
    yeast_gecko_biomass = "r_2111"
    yeast_gecko_succinate = "r_2056"
    gecko_eval_function = build_evaluation_function(
        "WYIELD",   # evaluating function id
        yeast_gecko_biomass,   # biomass reaction id
        yeast_gecko_succinate,   # target reaction id
        alpha = 0.3,                         # percentage of maximum target considered
        minBiomassValue = 0.03135  # minimum biomass for viable solution
        )

    # Final Step
    # Optimization
    gecko_optimization(
        model = yeast_single_pool,
        optim_type = "KO",   # KO or UO
        eval_func = gecko_eval_function,
        constraints = {},  # contraints format: {"reac_id": (lb, up), ...}
        critical_proteins = "auto",  # can also use a list of proteins
        isMultiProc = False,  # replace with True to use multi-processing
        size = 10,  # chose max number of alterations, None uses default
        output_file = "optimization_results.csv"   # chose path of the results file
        )
