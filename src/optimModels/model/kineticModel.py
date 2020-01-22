from collections import OrderedDict
from libsbml import readSBMLFromFile
from framed.io.sbml import _load_compartments, _load_metabolites, _load_reactions, _load_global_parameters, \
    _load_local_parameters, _load_ratelaws, _load_assignment_rules, _load_concentrations
from framed.model.odemodel import ODEModel
import re
from optimModels.utils.utils import MyTree


def load_kinetic_model(filename, kmap = None):
    """
    Load the dynamic model present in a SBML file.
    Args:
        filename (str): SBBML file.
        kmap (dict): Dictionary with the parameters that can be used in the strain
            optimization process for each reaction.{id_reaction: [param1, param2]}

    Returns: KineticModel
        Contains all information related with the dynamic model
            (reactions, kinetic equations, metabolites, compartments, etc.)
    """
    document = readSBMLFromFile(filename)
    sbmlModel = document.getModel()

    if sbmlModel is None:
        raise IOError('Failed to load model.')

    model = KineticModel(sbmlModel.getId())

    _load_compartments(sbmlModel, model)
    _load_metabolites(sbmlModel, model)
    _load_reactions(sbmlModel, model)
    _load_concentrations(sbmlModel, model)
    _load_global_parameters(sbmlModel, model)
    _load_local_parameters(sbmlModel, model)
    _load_ratelaws(sbmlModel, model)
    _load_assignment_rules(sbmlModel, model)

    # parse rates, rules and xdot expressions
    model._set_parsed_attr()

    if isinstance(kmap, str):
        aux = OrderedDict([(rId, re.findall(kmap, ratelaw)) for rId, ratelaw in model.ratelaws.items()])
        # aux = OrderedDict([(rId ,  [rId+"_"+x for x in re.findall("(rmax\w*)", ratelaw)])
        #                    for rId, ratelaw in model.ratelaws.items()])  # CHASSAGNOLE
        model.reacParamsFactors = OrderedDict([(rId, params) for rId, params in aux.items() if len(params) > 0])
    else:
        model.set_reactions_parameters_factors(kmap)
    return model


class KineticModel(ODEModel):
    """ Class to store information of dynamic models.
    This class is an extension of ODEModel class from FRAMED package.
        The methods  *build_ode* and *get_ode* were overrided to
    support the manipulations over the parameters or enzyme concentration level during the strain optimization process.
    The ode function returned by *get_ode* replace the negative concentrations,
        received as argument, by 0 if the concentration is smaller than
    absolute tolerance or raise an exception if the concentration is a significant negative value.
    """

    def __init__(self, modelId):
        """
        Create a instance of dynamicModel class.

        Args:
            modelId (str): A valid unique identifier
        """
        ODEModel.__init__(self, modelId)
        self.reacParamsFactors = None  # reaction->parameters association
        self.parsedRates = None
        self.parsedRules = None
        self.parsedXdot = None

    def _set_parsed_attr(self):
        self.parsedRates = {rId: self.parse_rate(rId, ratelaw) for rId, ratelaw in self.ratelaws.items()}
        aux = {pId: self.parse_rule(rule, self.parsedRates)
               for pId, rule in self.assignment_rules.items()}

        trees = [_build_tree_rules(vId, aux) for vId in aux.keys()]
        order = _get_oder_rules(trees)

        self.parsedRules = OrderedDict([(rule_id, aux[rule_id]) for rule_id in order])
        self.parsedXdot = {mId: self.print_balance(mId) for mId in self.metabolites}

    def build_ode(self, factors):
        """
        Args:
            factors (dict): The key is the parameter identifier and the value is the level of change values
                between 0 and 1 represent a under expression, above 1 a over
        expression and 0 to represent the knockouts.

        Returns: (str)  Returns  a string with the ode system.
        """

        # factors: ["vmax1": 0, "vmax2"=2, "ENZYME_ID":0]
        # divide vmax parameters from enzymes expression levels
        factorsEnz = OrderedDict([(k, v) for k, v in factors.items() if k in self.metabolites.keys()])
        factorsParam = OrderedDict([(k, v) for k, v in factors.items() if k not in factorsEnz.keys()])

        ruleExprs = ["    v['{}'] = {}".format(pId, self.parsedRules[pId])
                     for pId in self.parsedRules.keys()]

        rateExprs = []
        for rId in self.reactions.keys():
            newExp = self.parsedRates[rId]
            if rId in self.reacParamsFactors.keys():
                toModify = set(factorsParam.keys()).intersection(self.reacParamsFactors[rId])
                if len(toModify) > 0:
                    for elem in toModify:
                        newExp = re.sub(r"([pv]\['" + elem + r"'\])", str(factorsParam[elem]) + r" * \1", newExp)
            rateExprs.append("    r['{}'] = {}".format(rId, newExp))

        balances = []
        for m_id in self.metabolites.keys():
            exp = self.parsedXdot[m_id]
            if m_id in factorsEnz.keys():
                if factors[m_id] == 0:
                    newExp = "0"
                else:
                    newExp = re.sub(r"\+\s*(\d)", r"+ \1 * " + str(factorsEnz[m_id]), exp)
                balances.append(' ' * 8 + newExp)
            else:
                balances.append(' ' * 8 + exp)

        func_str = 'def ode_func(t, x, r, p, v):\n\n' + \
                   '\n'.join(ruleExprs) + '\n\n' + \
                   '\n'.join(rateExprs) + '\n\n' + \
                   '    dxdt = [\n' + \
                   ',\n'.join(balances) + '\n' + \
                   '    ]\n\n' + \
                   '    return dxdt\n'
        return func_str

    def get_ode(self, r_dict = None, params = None, factors = None):
        """
        Build the ODE system.

        Args:
            r_dict (dict): This variable is used to store the reaction rates.
            params (dict): Parameters and the new values used to replace
                the original parameters present in the SBML model.
            factors (dict): The key is the parameter identifier and the value is the level of change values between
            0 and 1 represent a under expression, above 1 a over expression and 0 to represent the knockouts.

        Returns:
            func: A function used to solve the ODE system.

        """

        p = self.merge_constants()
        v = self.variable_params.copy()

        if r_dict is not None:
            r = r_dict
        else:
            r = {}

        if params:
            p.update(params)

        exec(self.build_ode(factors), globals())
        ode_func = eval('ode_func')
        # print "Parameters"
        # print p
        # print v
        # print "-----"
        # print (self.build_ode(factors))

        # f = lambda t, x: ode_func(t, x, r, p, v)
        def f(t, x):
            return ode_func(t, x, r, p, v)
        return f

    def set_reactions_parameters_factors(self, kmap):
        """
        Set a new map with the parameters that can be changed for each reaction.

        Args:
            kmap (dict): The keys is the reaction identifier and the value a list of parameters
                which can be used to simulate modifications( KO, under/ over expression)

        """
        self.reacParamsFactors = OrderedDict(kmap)

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)


# auxiliary functions to set the assignment rules by the correct order in the ODE system
# TODO: Check if the framed in python 3.6 already correct the order of rules
def _build_tree_rules(parent, rules):
    regexp = r"v\[\'(.*?)\'\]"
    children = re.findall(regexp, rules[parent])
    if len(children) == 0:
        return MyTree(parent, None)
    else:
        childrenTrees = [_build_tree_rules(child, rules) for child in children]
        return MyTree(parent, childrenTrees)


def _get_oder_rules(trees):
    res = []
    for tree in trees:
        new_elems = _get_order_nodes(tree)
        [res.append(item) for item in new_elems if item not in res]
    return res


def _get_order_nodes(tree):
    res = [tree.name]
    if len(tree.children) > 0:
        for child in tree.children:
            res = _get_order_nodes(child) + res
    return res
