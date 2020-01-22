from optimModels.utils.constantes import solverStatus


class SimulationResult():
    """
    Represents the result of a metabolic model simulation at steady-state.
    """

    def __init__(self, model, solverStatus, ssFluxesDistrib,
                 overrideSimulProblem = None):
        """
        Create a Simulationresult instance.

        Args:
            modelId (str): Identification of metabolic model
            model: The model
            solverStatus (int): Simulation result (OPTIMAL = 0, UNKNOWN = 1, ERROR = 2).
            ssFluxesDistrib (dict): Fluxes distribution achieved in steady state.
            overrideSimulProblem (OverrideSimulProblem): Modifications over the metabolic model.
        """
        self.model = model
        self.solverStatus = solverStatus
        self.ssFluxesDistrib = ssFluxesDistrib
        self.overrideSimulProblem = overrideSimulProblem

    def get_solver_status(self):
        """
        Returns  the solver status result. (see optimModels.utils.constants.solverStatus)
        """
        return self.solverStatus

    def get_override_simul_problem(self):
        """
        Gets the override simulation problem.
        """
        return self.overrideSimulProblem

    def get_fluxes_distribution(self):
        """
        Gets the steady-state flux distribution {reactionId: fluxValue}.
        """
        return self.ssFluxesDistrib

    def get_model(self):
        """
        gets the steady-state model
        """
        return self.model


class StoicSimulationResult(SimulationResult):

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def print(self):
        print("Phenotype Simulation")
        print("------------------------")
        print("model id: " + self.model.id)
        print("status: " + solverStatus.get_status_str(self.solverStatus))
        print("fluxes: ")
        for k, v in self.ssFluxesDistrib.items():
            print("     " + k + " = " + str(v))

        if self.overrideSimulProblem:
            print("mofifications:")
            for k, v in self.overrideSimulProblem.get_modifications().items():
                print("     " + k + " = " + str(v))
        print("------------------------")

    # TODO: test
    def get_net_conversion(self, biomassId = None):
        '''
           Returs a string representation of the net conversion

           args:
                biosmassId (str) : optional
        '''

        ssFluxes = self.get_fluxes_distribution()

        left = ""
        right = ""

        # just for pretty print
        firstLeft, firstRight = True, True

        reactions = self.model.reactions
        for r_id in reactions.keys():
            fluxValue = ssFluxes[r_id]
            sub = reactions[r_id].get_substrates()
            prod = reactions[r_id].get_products()
            # if rId is a drain reaction
            if fluxValue != 0.0 and (len(sub) == 0 or len(prod) == 0):
                m = sub + prod
                if fluxValue < 0:
                    if firstLeft:
                        firstLeft = False
                    else:
                        left = left + " + "
                    left = left + str(-1 * fluxValue)
                    left = left + " " + m[0]
                else:
                    if firstRight:
                        firstRight = False
                    else:
                        right = right + " + "
                    right = right + str(fluxValue)
                    right = right + " " + m[0]

        if biomassId and biomassId in ssFluxes.keys():
            biomassFlux = ssFluxes[biomassId]
            if biomassFlux > 0:
                if firstRight:
                    firstRight = False
                else:
                    right = right + " + "
                right = right + str(biomassFlux)
                right = right + " " + biomassId

        return left + " --> " + right


class GeckoSimulationResult(SimulationResult):

    def __init__(self, model, solverStatus, ssFluxesDistrib = None, protConcentrations = None,
                 overrideSimulProblem = None):
        super().__init__(model, solverStatus, ssFluxesDistrib, overrideSimulProblem)
        self.protConcentrations = protConcentrations

    def get_protein_concentrations(self):
        """
        Gets the protein concentrations in steady-state {proteinId: concentration value}.
        """
        return self.protConcentrations

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def print(self):
        print("Phenotype Simulation")
        print("------------------------")
        print("model id: " + self.model.id)
        print("status: " + solverStatus.get_status_str(self.solverStatus))
        print("fluxes: ")
        for k, v in self.ssFluxesDistrib.items():
            print("     " + k + " = " + str(v))

        print("protein concentrations: ")
        for k, v in self.protConcentrations.items():
            print("     " + k + " = " + str(v))

        if self.overrideSimulProblem:
            print("mofifications:")
            for k, v in self.overrideSimulProblem.get_modifications().items():
                print("     " + k + " = " + str(v))
        print("------------------------")


class kineticSimulationResult(SimulationResult):
    """ Represents the result of a dynamic metabolic model simulation on steady-state.

    Args:
        modelId (str): identification of metabolic model
        solverStatus (int): simulation result (OPTIMAL = 0, UNKNOWN = 1, ERROR = 2).
        ssFluxesDistrib (dict): fluxes distribution achieved in steady state.
        ssConcentrations (dict): metabolites concentration in steady state.
        overrideSimulProblem (overrideKineticSimulProblem): modifications over the metabolic model.
    """

    def __init__(self, model, solverStatus, ssFluxesDistrib, ssConcentrations = None,
                 overrideSimulProblem = None):
        self.ssConcentrations = ssConcentrations
        super().__init__(model, solverStatus, ssFluxesDistrib, overrideSimulProblem)

    def get_steady_state_concentrations(self):
        """
        Gets the metabolite concentrations in steady-state {metaboliteId: concentration value}.
        """
        return self.ssConcentrations

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def print(self):
        print("Phenotype Simulation")
        print("------------------------")
        print("model id: " + self.model.id)
        print("status: " + solverStatus.get_status_str(self.solverStatus))
        print("fluxes: ")
        for k, v in self.ssFluxesDistrib.items():
            print("     " + k + " = " + str(v))

        print("concentrations: ")
        for k, v in self.ssConcentrations.items():
            print("     " + k + " = " + str(v))

        if self.overrideSimulProblem:
            print("mofifications:")
            for k, v in self.overrideSimulProblem.get_modifications().items():
                print("     " + k + " = " + str(v))
        print("------------------------")
