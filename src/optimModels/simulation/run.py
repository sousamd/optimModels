from optimModels.simulation.simul_problems import KineticSimulationProblem
from optimModels.simulation.override_simul_problem import OverrideKineticSimulProblem


def kinetic_simulation(model, parameters = None, factors = None, time = 1e9):
    """
    Runs a phenotype simulation using dynamic models.
    Args:
        model (kineticModel): The kinetic metabolic model.
        parameters (dict): List of parameters that will be set with new values (ex: Dilution, initial concentrations).
        factors (dict): Values to by multiplied to the vMax parameters
            (KO: the value should be 0, Under: value between 0 and 1, Over: value higher than 1)
        time (float): End time for steady-state.

    Returns (kineticSimulationResults): The function returns the best solutions found in strain optimization.
        The kineticSimulationResults have the
        flux distribution and metabolites concentration on steady-state, and the modifications made over the
        original model.
    """

    simulProblem = KineticSimulationProblem(model, parameters = parameters, tSteps=[0, time], timeout=None)

    if factors:
        override = OverrideKineticSimulProblem(factors = factors)
        result = simulProblem.simulate(override)
    else:
        result = simulProblem.simulate()

    return result
