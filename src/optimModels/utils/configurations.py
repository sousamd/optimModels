from optimModels.utils.constantes import solverMethod


class SolverConfigurations:
    ABSOLUTE_TOL = 1e-9
    RELATIVE_TOL = 1e-6
    N_STEPS =10000

class StoicConfigurations:
    SOLVER = 'cplex'
    SOLVER_METHOD = 'FBA'
    DEFAULT_LB = -99999
    DEFAULT_UB = 99999

class GeckoConfigurations:
    SCALE_CONSTANT = 1000000
    SOLVER_METHOD = 'FBA'

class KineticConfigurations:
    STEADY_STATE_TIME = 1e20
    SOLVER = "odespy"
    SOLVER_METHOD = solverMethod.LSODA  # ode solver method used in the phenotype simulation
    SOLVER_TIMEOUT = 6000  # maximum time allowed by simulation




