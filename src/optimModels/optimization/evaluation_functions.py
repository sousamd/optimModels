from abc import ABCMeta, abstractmethod
from itertools import chain
from optimModels.utils.configurations import StoicConfigurations
import optimModels.simulation.analysis as an
from framed.cobra.simulation import FBA, pFBA, MOMA, lMOMA, ROOM
from framed.solvers import solver_instance


class EvaluationFunction:
    """
    This abstract class should be extended by all evaluation functions classes.

    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_fitness(self, simulationResult, candidate):
        return

    @abstractmethod
    def method_str(self):
        return

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)


class MinCandSizeAndMaxTarget(EvaluationFunction):
    """
    This evaluation function finds the solution with the minimum candidate size and levels considering the maximization
    of the target flux.

    Args:
        maxCandidateSize (int): maximum of candidate size
        maxTargetFlux (str): reaction id to maximize

    """

    def __init__(self, maxCandidateSize, maxTargetFlux):
        self.maxCandidateSize = maxCandidateSize
        self.objective = maxTargetFlux

    def get_fitness(self, simulResult, candidate):
        fluxes = simulResult.get_fluxes_distribution()
        numModifications = len(list(chain.from_iterable(candidate)))
        sumObj = 0
        for rId, ub in self.objective.items():
            ub = StoicConfigurations.DEFAULT_UB if ub is None else ub
            f = 1 if fluxes[rId] >= ub else 1 - ((ub - fluxes[rId]) / ub)
            sumObj = sumObj + f
        objFactor = sumObj / len(self.objective)

        return objFactor / numModifications

    def method_str(self):
        return "Minimize the number of modifications while maximize the target flux."

    @staticmethod
    def get_id():
        return "MinNumberReacAndMaxFlux"

    @staticmethod
    def get_name():
        return "Minimize the number of modifications while maximize the target flux."

    @staticmethod
    def get_parameters_ids():
        return ["Maximum of modifications allowed", "Target reactions"]


class MinCandSizeWithLevelsAndMaxTarget(EvaluationFunction):
    # TODO: Validar o que acontece em caso do ub do target ser 0 ou seja o fluxes[rId] é negativo (ver se há hipotese de isto acontecer)

    def __init__(self, maxCandidateSize, levels, maxTargetFlux):
        self.maxCandidateSize = maxCandidateSize
        self.levels = levels
        self.objective = maxTargetFlux

    def get_fitness(self, simulResult, candidate):
        fluxes = simulResult.get_fluxes_distribution()
        maxUptake = len(candidate) * self.levels[-1]
        sumUptake = 0
        sumObj = 0

        for rId in candidate.keys():
            sumUptake = sumUptake + candidate[rId]

        for rId, ub in self.objective.items():
            ub = StoicConfigurations.DEFAULT_UB if ub is None else ub
            f = 1 if fluxes[rId] >= ub else 1 - ((ub - fluxes[rId]) / ub)
            sumObj = sumObj + f
        objFactor = sumObj / len(self.objective)

        # best solution when factors are close to 1
        upFactor = sumUptake / maxUptake
        lenFactor = len(candidate) / self.maxCandidateSize

        return objFactor / (upFactor * lenFactor)

    def method_str(self):
        return "Minimize the number and the fluxes of candidate while maximize the target flux."

    @staticmethod
    def get_id():
        return "MinCandSizeWithLevelsAndMaxTarget"

    @staticmethod
    def get_name():
        return "Minimize the number and the fluxes of candidate while maximize the target flux."

    @staticmethod
    def get_parameters_ids():
        return ["Maximum number of modifications allowed", "Levels", "Target reactions"]


class MinCandSize(EvaluationFunction):
    """
    This class implements the "minimization of number of reactions" objective function. The fitness is given by
    1 - size(candidate)/ max_candidate_size, where the max_candidate_size is the maximum size that a candidate can have
    during optimization.

    Args:
        maxCandidateSize(int): Maximum size allowed for candidate
        minFluxes (dict): Minimal value for fluxes to consider fitness different of 0 (key: reaction id, value: minimum of flux).

    """

    def __init__(self, maxCandidateSize, minFluxes):
        self.maxCandidateSize = maxCandidateSize
        self.minFluxes = minFluxes

    def get_fitness(self, simulResult, candidate):
        fluxes = simulResult.get_fluxes_distribution()
        for rId, flux in self.minFluxes.items():
            if fluxes[rId] < flux:
                return 0
        return 1 - len(candidate) / (self.maxCandidateSize + 1)

    def method_str(self):
        return "Minimum number of active reactions."

    @staticmethod
    def get_id():
        return "MinNumberReac"

    @staticmethod
    def get_name():
        return "Minimum number of active reactions."

    @staticmethod
    def get_parameters_ids():
        return ["Number maximum of modification allowed", "Minimum of targets flux values."]


class TargetFlux(EvaluationFunction):
    """
    This class implements the "target flux" objective function. The fitness is given by the flux value of the target reaction.

    Args:
        targetReactionId (str): Reaction identifier of the target compound production.

    """

    def __init__(self, targetReactionId):
        # TODO: take only the first element - done
        self.targetReactionId = targetReactionId[0]

    def get_fitness(self, simulResult, candidate):
        fluxes = simulResult.get_fluxes_distribution()
        if self.targetReactionId not in list(fluxes.keys()):
            raise ValueError("Reaction id is not present in the fluxes distribution.")
        return fluxes[self.targetReactionId]

    def method_str(self):
        return "Target Flux: " + self.targetReactionId

    @staticmethod
    def get_id():
        return "targetFlux"

    @staticmethod
    def get_name():
        return "Target Flux"

    @staticmethod
    def get_parameters_ids():
        return ["Target reaction id"]


class BPCY(EvaluationFunction):
    """
    This class implements the "Biomass-Product Coupled Yield" objective function. The fitness is given by the equation:
    (biomass_flux * product_flux)/ uptake_flux

    Args:
        biomassId (str): Biomass reaction identifier
        productId (str): Target product reaction identifier
        uptakeId (str): Reaction of uptake

    If no substract is defined, ie uptakeId is None, a substracted value of 1.0 is considered
        """

    def __init__(self, biomassId, productId, uptakeId, **kwargs):
        self.biomassId = biomassId
        self.productId = productId
        self.uptakeId = uptakeId
        self.method = kwargs.get('method', None)

    def get_fitness(self, simulResult, candidate):
        if self.method is not None:
            model = simulResult.get_model()
            _constraints = {}
            simulResult.get_override_simul_problem().get_modifications()
            if self.method == 'FBA':
                solution = FBA(model, constraints = _constraints)
            elif self.method == 'pFBA':
                solution = pFBA(model, constraints = _constraints)
            elif self.method == 'MOMA':
                solution = MOMA(model, constraints = _constraints)
            elif self.method == 'lMOMA':
                solution = lMOMA(model, constraints = _constraints)
            elif self.method == 'ROOM':
                solution = ROOM(model, constraints = _constraints)
            else:
                raise ValueError("Unknown method ", self.method)
            ssFluxes = solution
        else:
            ssFluxes = simulResult.get_fluxes_distribution()

        ids = list(ssFluxes.keys())
        if self.biomassId not in ids or self.productId not in ids or self.uptakeId not in ids:
            raise ValueError(
                "Reaction ids is not present in the fluxes distribution. Please check id objective function is correct.")
        if abs(ssFluxes[self.uptakeId]) == 0:
            return 0
        else:
            return (ssFluxes[self.biomassId] * ssFluxes[self.productId]) / abs(ssFluxes[self.uptakeId])

    def method_str(self):
        return "BPCY =  (" + self.biomassId + " * " + self.productId + ") / " + self.uptakeId

    @staticmethod
    def get_id():
        return "BPCY"

    @staticmethod
    def get_name():
        return "Biomass-Product Coupled Yield"

    @staticmethod
    def get_parameters_ids():
        return ["Biomass id", "Product id", "Uptake id"]


class BP_MinModifications(EvaluationFunction):
    """
        This class is based the "Biomass-Product Coupled Yield" objective function but considering the candidate size. The fitness is given by the equation:
        (biomass_flux * product_flux)/ candidate_size)

        Args:
            biomassId (str): biomass reaction identifier
            productId (str): target product reaction identifier

        """

    def __init__(self, biomassId, productId):
        self.biomassId = biomassId
        self.productId = productId

    def get_fitness(self, simulResult, candidate):
        ssFluxes = simulResult.get_fluxes_distribution()
        ids = list(ssFluxes.keys())
        if self.biomassId not in ids or self.productId not in ids:
            raise ValueError(
                "Reaction ids is not present in the fluxes distribution. Please check id objective function is correct.")
        size = len(candidate)
        # optimization of medium and KO .. the candidate is a tuple of lists
        if (isinstance(candidate[0], list)):
            size = len(candidate[0]) + len(candidate[1])
        print(candidate, str(ssFluxes[self.biomassId]), str(ssFluxes[self.productId]))
        return (ssFluxes[self.biomassId] * ssFluxes[self.productId]) / size

    def method_str(self):
        return "BP_MinModifications=  (" + self.biomassId + " * " + self.productId + ") / candidate_size "

    @staticmethod
    def get_id():
        return "BP_MinModifications"

    @staticmethod
    def get_name():
        return "Biomass-Product with minimun of modifications"

    @staticmethod
    def get_parameters_ids():
        return ["Biomass id", "Product id"]


class WYIELD(EvaluationFunction):
    """
    Weighted Yield (WYIELD) objective function

    Arguments:
            biomassId (str): biomass reaction identifier
            productId (str): target product reaction identifier
    **kwargs options:
            minBiomassValue (float): minimum biomass value (default None, in which case the minBiomassPercentage is used)
            minBiomassPercentage (float) : Instead of defining explicitly a minimum biomass value, a percentage of the
                                           wild type biomass is used. Only used when no minBiomassValue is defined
                                           (default minBiomassPercentage:0.10)
            alpha (float): tradeoff between the min and max FVA of the target product.
                           Must be in range [0,1]  (default alpha:0.5)
            scale (boolean): defines if the WYIELD is devided by the biomass of the simulated result (default false)
    """

    def __init__(self, biomassId, productId, **kwargs):
        self.biomassId = biomassId
        self.productId = productId

        # kwargs parameters
        self.minBiomassValue = kwargs.get('minBiomassValue', None)
        self.minBiomassPercentage = kwargs.get('minBiomassPercentage', 0.10)
        self.alpha = kwargs.get('alpha', 0.5)
        self.scale = kwargs.get('scale', False)

        if self.alpha > 1 or self.alpha < 0:
            raise ValueError("The value of the tradeoff parameter alpha should be in range 0 to 1.")

        self.worstFitness = 0.0

    def get_fitness(self, simulResult, candidate):
        fvaMaxProd = 0.0
        fvaMinProd = 0.0

        model = simulResult.get_model()

        ssFluxes = simulResult.get_fluxes_distribution()
        ids = list(ssFluxes.keys())
        if self.biomassId not in ids or self.productId not in ids:
            raise ValueError(
                "Reaction ids are not present in the fluxes distribution. Please check if the objective function ids are correct.")

        biomassFluxValue = ssFluxes[self.biomassId] * 0.99999

        constraints = simulResult.get_override_simul_problem().get_modifications()

        # todo: necessario calcular isto de novo? ASK
        if self.minBiomassValue is None or self.minBiomassValue < 0.0:
            solution = an.FBA(model, objective = {self.biomassId: 1}, constraints = constraints)
            wtBiomassValue = solution[self.biomassId]

            minBiomass = wtBiomassValue * self.minBiomassPercentage
            self.minBiomassValue = minBiomass
        else:
            minBiomass = self.minBiomassValue

        # add biomass constraint
        constraints[self.biomassId] = (biomassFluxValue, 100000.0)

        # only need to simulate FVA min if alpha is lesser than 1, otherwise it will always be zero
        if (self.alpha < 1):
            fvaMinResult = an.FBA(model, objective = {self.productId: 1}, minimize = True, constraints = constraints)
            if fvaMinResult:
                fvaMinProd = fvaMinResult[self.productId]
                print('min', fvaMinProd)

        # only need to simulate FVA max if alpha is larger than 0, otherwise it will always be zero
        if (self.alpha > 0):
            fvaMaxResult = an.FBA(model, objective = {self.productId: 1}, constraints = constraints)
            if fvaMaxResult:
                fvaMaxProd = fvaMaxResult[self.productId]
                print('max', fvaMaxProd)

        res = self.worstFitness

        if biomassFluxValue > minBiomass:
            res = (self.alpha * fvaMaxProd + (1 - self.alpha) * fvaMinProd)
            if self.scale:
                res = res / biomassFluxValue
        return res, fvaMaxProd, fvaMinProd

    def method_str(self):
        return "WYIELD "

    @staticmethod
    def get_id():
        return "WYIELD"

    @staticmethod
    def get_name():
        return "Weighted Yield"

    @staticmethod
    def get_parameters_ids():
        return ["Biomass id", "Product id"]


def build_evaluation_function(id, *args, **kwargs):
    """
    Function to return an evaluation function instance.

    Args:
        id (str): Name of the objective function. The implemented objective functions should be registed in constants.objFunctions class
        *args (list of str): The number of arguments depends of the objective function chosen by user.
        **kwargs dictionary of additional optional arguments
    Returns:
        EvaluationFunction: return an evaluation function instance.
    """

    if id == BPCY.get_id():
        objFunc = BPCY(args[0], args[1], args[2])
    elif id == WYIELD.get_id():
        objFunc = WYIELD(args[0], args[1], **kwargs)
    elif id == TargetFlux.get_id():
        objFunc = TargetFlux(args[0])
    elif id == MinCandSize.get_id():
        objFunc = MinCandSize(args[0], args[1])
    elif id == BP_MinModifications.get_id():
        objFunc = BP_MinModifications(args[0], args[1])
    elif id == MinCandSizeWithLevelsAndMaxTarget.get_id():
        objFunc = MinCandSizeWithLevelsAndMaxTarget(args[0], args[1], args[2])
    elif id == MinCandSizeAndMaxTarget.get_id():
        objFunc = MinCandSizeAndMaxTarget(args[0], args[1])
    else:
        raise Exception("Wrong objective function!")

    return objFunc

