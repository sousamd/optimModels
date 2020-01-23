from framed.io.sbml import load_cbmodel
from optimModels.simulation.simul_problems import StoicSimulationProblem
from cobra.io import read_sbml_model
from cobra.core import Reaction
from optimModels.optimization.decoders import DecoderReacUnderOverExpression
from optimModels.unittests.cbm_optim import get_nontargets

##################################################################################
# File used to validate some of the results obtained in the optimization process #
##################################################################################

# # iMM904
# # overexpression
# pyruvate_carboxilase = ["R_PC"] # cyto
# phophoenolpyruvate_carboxykinase = ["R_PPCK"] # cyto
# malate_dehydrogenase = ["R_MDH"]  # cyto         "R_MDHm" # mito        "R_MDHp" # pero
# fumarase = ["R_FUM"] # cyto  "R_FUMm" # mito
# fumarate_reductase = ["R_FRDcm"] # cyto(fum->succ)/mito(nadh)       "R_FRDm"  # mito
# isocitrate_lyase = ["R_ICL"] # cyto
# malate_synthase = ["R_MALS"] # cyto      "R_MALSp" # pero
#
# # knockout
# pyruvate_decarboxylase = ["R_PYRDC"] # cyto      "R_PYRDC2"  #cyto
# alcohol_dehydrogenase = ["R_ALCD2ir"] # cyto     "R_ALCD2irm"  # mito
#
# # transporters; unneeded transporters, everything in cytosol
# dicarboxylic_acid_transporter = {
#     "alpha-ketoglutarate": ["R_AKGMAL"],
#     "aspartate": ["R_ASPt7", "R_CBASPtn", "R_ASPGLU2m", "R_ASPt2m", "R_ASPt2n", "R_ASPt2r", "R_ASPt5n"],
#     "malate": ["R_MALtm", "R_MALt2r", "R_3C3HMPt", "R_3C3HMPtm", "R_AKGMAL"],
#     "fumarate": ["R_SUCFUMtm", "R_FUMt2r"],
#     "oxaloacetate": ["R_OAAt", "R_OAAt2m"],
#     "succinate": ["R_SUCCt2r", "R_SUCCtm", "R_SUCFUMtm"],
#     "glutamate": ["R_ASPGLU2m", "R_GLUt2r", "R_GLUt2n", "R_GLUt5m",
#                   "R_GLUt7", "R_GLUt7m", "R_E4HGLUtm", "R_E4HGLUtp"],
#     "glutarate": ["R_4H2OGLTtm", "R_4H2OGLTtp", "R_AKGt2n", "R_AKGt2r"]
#     }
# over_prots = [pyruvate_carboxilase, phophoenolpyruvate_carboxykinase, malate_synthase,
#               malate_dehydrogenase, fumarate_reductase, fumarase, isocitrate_lyase]
# ko_prots = [pyruvate_decarboxylase, alcohol_dehydrogenase]


# yeastGEM
# overexp
pyruvate_carboxylase = ["r_0958"]
phosphoenolpyruvate_carboxykinase = ["r_0884"]
malate_dehydrogenase = ["r_0714"]  # , "r_0713"]
fumarase = ["r_0452"]  # , "r_0451"]
fumarate_reductase = ["r_1000"]  # "r_0454", "r_0455",
isocitrate_lyase = ["r_0662"]
malate_synthase = ["r_0716"]  # , "r_0717"]

# ko
alcohol_dehydrogenase = ["r_0163", "r_0165", "r_2115"]
pyruvate_decarboxylase = ["r_0959"]
koprots = alcohol_dehydrogenase + pyruvate_decarboxylase

# extra6
citrate_synthase = ["r_0300"]  # , "r_0301"]
cit_to_cisaco = ["r_0302"]
cisaco_to_isocit = ["r_0280"]
isocitrate_dehydrogensase = ["r_0658", "r_2131"]


# Simulation:
def simul_iMM904SL(KO, constraints = None, withCobraPy = False):
    SBML_FILE = '../../../examples/models/iMM904SL_v6.xml'
    model = load_cbmodel(SBML_FILE, flavor = 'cobra')
    model.set_reaction_objective('R_biomass_SC5_notrace', coeff = 1)

    criticalReacs = get_nontargets('../../../examples/models/nontargets_iMM904SL_v6_UO.txt')

    simulProblem = StoicSimulationProblem(model, constraints = constraints)

    # LEVELS = [1e-3, 1e-2, 1e-1, 0.5, 0, 1, 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
    LEVELS = [0, 2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5]

    reactions = [x for x in simulProblem.get_internal_reactions() if
                 x not in criticalReacs and x not in simulProblem.objective.keys()]
    reactions.sort()

    from optimModels.optimization.decoders import DecoderReacUnderOverExpression, DecoderReacKnockouts
    decoder = DecoderReacUnderOverExpression(reactions, levels = LEVELS)
    # #decoder = DecoderReacKnockouts(reactions)

    candidate = decoder.decode_candidate_ids_to_index(identifiers = KO)
    override = decoder.get_override_simul_problem(candidate = candidate, simulProblem = simulProblem)

    res = simulProblem.simulate(override)
    res.print()

    print(res.ssFluxesDistrib['R_EX_succ_e'])
    print(res.ssFluxesDistrib['R_biomass_SC5_notrace'])

    # print((res.ssFluxesDistrib['R_EX_succ_e'] * res.ssFluxesDistrib['R_biomass_SC5_notrace']) / (
    #       -res.ssFluxesDistrib['R_EX_glc_e']))


def simul_iMM904(UO, constraints = None, withCobraPy = False):
    SBML_FILE = '../../../examples/models/iMM904.xml'
    model = load_cbmodel(SBML_FILE, flavor = 'bigg', load_gprs = None)
    # simulProblem = StoicSimulationProblem(model, constraints = constraints, objective = {"R_EX_succ_e": 1})
    simulProblem = StoicSimulationProblem(model, constraints = constraints, objective = {"R_BIOMASS_SC5_notrace": 1})

    criticalReacs = []

    LEVELS = [0, 2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5]
    reactions = [x for x in simulProblem.get_internal_reactions() if
                 x not in criticalReacs and x not in simulProblem.objective.keys()]
    reactions.sort()

    metabolite = "M_o2_r"
    print("Consumed: ", simulProblem.model.get_metabolite_consumers(metabolite))

    if not UO:
        res = simulProblem.simulate()
        res.print()
        print(res.ssFluxesDistrib["R_EX_succ_e"])
        print(res.ssFluxesDistrib["R_BIOMASS_SC5_notrace"])
        return res

    decoder = DecoderReacUnderOverExpression(reactions, levels = LEVELS)
    candidate = decoder.decode_candidate_ids_to_index(identifiers = UO)
    override = decoder.get_override_simul_problem(candidate = candidate, simulProblem = simulProblem)

    res = simulProblem.simulate(overrideSimulProblem = override)

    res.print()
    print(res.ssFluxesDistrib["R_EX_succ_e"])
    print(res.ssFluxesDistrib["R_BIOMASS_SC5_notrace"])

    return res


def simul_yeastGEM(UO, constraints = None, aeroflag = True, cobraflag = False):
    SBML_FILE = '../../../examples/models/yeastGEMv8_46_4.xml'
    biomass = "r_4041"
    succ = "r_2056"

    model = (load_cbmodel(SBML_FILE, flavor = 'fbc2'), read_sbml_model(SBML_FILE))[cobraflag]

    if not aeroflag:
        if not cobraflag:  # framed sรณ funciona anaerobiosis com a reacao de cima, cobrapy coma as duas
            model.add_reaction_from_str(
                "r_4598_alt: 0.000190000006114133 s_0529__91__c__93__ + "
                "9.99999974737875e-06 s_0687__91__c__93__ + "
                "0.000989999971352518 s_1405__91__c__93__ + "
                "1.20000004244503e-06 s_1475__91__c__93__ + "
                "6.34000025456771e-05 s_1487__91__c__93__ --> s_4205__91__c__93__")
            # # gecko appendix
            # model.add_reaction_from_str("r_4598_alt: 0.000190000006114133 s_0529__91__c__93__ + "
            #                             "9.99999974737875e-06 s_0687__91__c__93__ + "
            #                             "0.00264999992214143 s_1198__91__c__93__ + "
            #                             "0.000150000007124618 s_1203__91__c__93__ + "
            #                             "0.000569999974686652 s_1207__91__c__93__ + "
            #                             "0.00270000007003546 s_1212__91__c__93__ + "
            #                             "0.000989999971352518 s_1405__91__c__93__ + "
            #                             "1.20000004244503e-06 s_1475__91__c__93__ + "
            #                             "6.34000025456771e-05 s_1487__91__c__93__ --> s_4205__91__c__93__")
        else:
            model.add_reactions(
                [Reaction(id = "r_4598_alt", name = "r_4598_alt", subsystem = "c", upper_bound = 1000)]
                )
            model.reactions.r_4598_alt.add_metabolites(
                {"s_0529[c]": -0.000190000006114133, "s_0687[c]": -9.99999974737875e-06,
                 "s_1405[c]": -0.000989999971352518, "s_1475[c]": -1.20000004244503e-06,
                 "s_1487[c]": -6.34000025456771e-05, "s_4205[c]": 1.0})
            # gecko appendix
            model.reactions.r_4598_alt.add_metabolites(
                {"s_0529[c]": -0.000190000006114133, "s_0687[c]": -9.99999974737875e-06,
                 "s_1198[c]": -0.00264999992214143, "s_1207[c]": -0.000569999974686652,
                 "s_1203[c]": -0.000150000007124618, "s_1212[c]": -0.00270000007003546,
                 "s_1405[c]": -0.000989999971352518, "s_1475[c]": -1.20000004244503e-06,
                 "s_1487[c]": -6.34000025456771e-05, "s_4205[c]": 1.0})

    simulProblem = StoicSimulationProblem(model, constraints = constraints,
                                          objective = {biomass: 1}, method = "pFBA", minimize = False,
                                          withCobraPy = cobraflag)
    # simulProblem = StoicSimulationProblem(model, constraints = constraints,
    #                                       objective = {succ: 1}, method = "pFBA", minimize = False,
    #                                       withCobraPy = cobraflag)

    criticalReacs = []

    LEVELS = [0, 2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5]
    reactions = [x for x in simulProblem.get_internal_reactions() if
                 x not in criticalReacs and x not in simulProblem.objective.keys()]
    reactions.sort()

    if not UO:
        res = simulProblem.simulate()
    else:
        decoder = DecoderReacUnderOverExpression(reactions, levels = LEVELS)
        candidate = decoder.decode_candidate_ids_to_index(identifiers = UO)
        override = decoder.get_override_simul_problem(candidate = candidate, simulProblem = simulProblem)
        res = simulProblem.simulate(overrideSimulProblem = override)

    # res.print()
    print("succ: ", res.ssFluxesDistrib[succ])
    print("biom: ", res.ssFluxesDistrib[biomass])
    return res


if __name__ == '__main__':
    # #iMM904
    # UO = []
    # for prot in over_prots:
    #     UO += prot
    # KO = []
    # for prot in ko_prots:
    #     KO += prot
    #
    # UO_dic = {x: 2 for x in UO}
    #
    # constraints = {"R_BIOMASS_SC5_notrace": (0.026, 999999), "R_EX_o2_e": (-0.0015, 999999)}
    # constraints.update({x: (0, 0) for x in pyruvate_decarboxylase})
    #
    # res = simul_iMM904(UO = None, constraints = constraints)
    # # print("R_PC: ", res.ssFluxesDistrib["R_PC"])
    # # print("R_PPCK: ", res.ssFluxesDistrib["R_PPCK"])
    # # print("R_MDH: ", res.ssFluxesDistrib["R_MDH"])
    # # print("R_FUM: ", res.ssFluxesDistrib["R_FUM"])
    # # print("R_FRDcm: ", res.ssFluxesDistrib["R_FRDcm"])
    # # print("R_ICL: ", res.ssFluxesDistrib["R_ICL"])
    # # print("R_MALS: ", res.ssFluxesDistrib["R_MALS"])
    # # print("--------")
    # # print("R_EX_succ_e: ", res.ssFluxesDistrib["R_EX_succ_e"])
    #
    # Consumed = ['R_34HPPOR', 'R_3HAO', 'R_3OPHB5Hm', 'R_AACTOOR', 'R_APRTO2', 'R_ARAB14LO', 'R_C22STDS', 'R_C22STDSx',
    #            'R_C4STMO1', 'R_C4STMO2', 'R_C5STDS', 'R_CERH124_copy1', 'R_CERH126_copy1', 'R_CERH124_copy2',
    #            'R_CERH126_copy2', 'R_CERS324', 'R_CERS326', 'R_CHLSTI', 'R_CPPPGO', 'R_DESAT14', 'R_DESAT16',
    #            'R_DESAT18', 'R_DESAT18_2', 'R_DHORDi', 'R_FAS141', 'R_FAS161', 'R_FAS181', 'R_HGNTOR', 'R_KYN3OX',
    #            'R_LNS14DM', 'R_LNS14DMx', 'R_O2ter', 'R_O2tm', 'R_PDX5POi', 'R_POLYAO', 'R_POLYAO2', 'R_POLYAO3',
    #            'R_PSPHS', 'R_PYAM5PO', 'R_PYDXNO', 'R_PYDXO', 'R_TAUDO', 'R_TRPO2']
    #
    # for reac in Consumed:
    #     if res.ssFluxesDistrib[reac]:
    #         print(reac, res.ssFluxesDistrib[reac])

    # yeastGEM
    aeroflag = True
    # aeroflag = False
    aero_biom = 0.9032244521605517  # 0.9064423393108797
    anaero_biom = 0.19552584273269918  # 0.19788792994974622
    biom = (anaero_biom, aero_biom)[aeroflag]
    biom_percent = 0.50

    constraints = {
        "r_4041": (biom * biom_percent, 1000),  # biomass
        "r_1714": (-10, 1000),  # glucose
        **{x[0]: (0, 1000) for x in [
            pyruvate_carboxylase, phosphoenolpyruvate_carboxykinase,
            fumarate_reductase, isocitrate_lyase, malate_synthase]},  # force flux in direction
        **{x[0]: (-1000, -0) for x in [malate_dehydrogenase, fumarase]},  # force flux in direction
        }

    constraints_anaero = {
        # # gecko appendix
        # "r_1915": (-1000, 1000),  # lanosterol
        # "r_2134": (-1000, 1000),  # demethyllanosterol
        # "r_2137": (-1000, 1000),  # ergosta-tetraen-beta-ol
        # "r_2106": (-1000, 1000),  # zymosterolnsase
        # "r_0472": (-1000, 0),  # glutamate synthase
        # "r_4046": (0, 0),  # NGAM

        "r_1992": (0, 1000),  # oxygen
        "r_1994": (-1000, 1000),  # palmitoleate
        "r_2189": (-1000, 1000),  # oleate
        "r_1757": (-1000, 1000)  # ergosterol
        }

    if not aeroflag:
        constraints.update(constraints_anaero)

    # Overexps
    overprots = pyruvate_carboxylase + phosphoenolpyruvate_carboxykinase\
        + malate_dehydrogenase + fumarase + fumarate_reductase + isocitrate_lyase + malate_synthase

    UO_dic = {x: 2 for x in overprots}

    # KOs
    # constraints.update({prot: (0, 0) for prot in alcohol_dehydrogenase})
    # constraints.update({prot: (0, 0) for prot in pyruvate_decarboxylase})

    res = simul_yeastGEM(UO = None, constraints = constraints, aeroflag = aeroflag, cobraflag = False)

    if res:
        print("-~-~-~-~-~-")
        print("pyruvate carboxylase:", res.ssFluxesDistrib[pyruvate_carboxylase[0]])
        print("phosphoenolpyruvate carboxykinase:", res.ssFluxesDistrib[phosphoenolpyruvate_carboxykinase[0]])
        print("malate dehydrogenase:", res.ssFluxesDistrib[malate_dehydrogenase[0]])
        print("fumarase:", res.ssFluxesDistrib[fumarase[0]])
        print("fumarate reductase:", res.ssFluxesDistrib[fumarate_reductase[0]])
        print("isocitrate lyase:", res.ssFluxesDistrib[isocitrate_lyase[0]])
        print("malate synthase:", res.ssFluxesDistrib[malate_synthase[0]])
        print("oxygen:", res.ssFluxesDistrib["r_1714"])

        # print("alcohol dehydrogenase:",  res.ssFluxesDistrib[alcohol_dehydrogenase[0]])
        # print("alcohol dehydrogenase:",  res.ssFluxesDistrib[alcohol_dehydrogenase[1]])
        # print("alcohol dehydrogenase:",  res.ssFluxesDistrib[alcohol_dehydrogenase[2]])
        # print("pyruvate decarboxylase:",  res.ssFluxesDistrib[pyruvate_decarboxylase[0]])
        # # extra6
        # print("citrate synthase:",  res.ssFluxesDistrib[citrate_synthase[0]])
        # print("citrate to cis-aconitate:",  res.ssFluxesDistrib[cit_to_cisaco[0]])
        # print("cis-aconitate to isocitrate:",  res.ssFluxesDistrib[cisaco_to_isocit[0]])
        # print("isocitrate dehydrogenase:",  res.ssFluxesDistrib[isocitrate_dehydrogensase[0]])
        # print("isocitrate dehydrogenase:",  res.ssFluxesDistrib[isocitrate_dehydrogensase[1]])
        # # print("4-aminobutyrate--2-oxoglutarate transaminase:",  res.ssFluxesDistrib[pyruvate_decarboxylase[0]])
        #
        # # where did the citrate go
        # print("AKG tranporter: ", res.ssFluxesDistrib["r_1112"])
        # print("citrate transport: ", res.ssFluxesDistrib["r_1126"])
        # print("citrate transport: ", res.ssFluxesDistrib["r_1127"])
        # print("citrate transport: ", res.ssFluxesDistrib["r_1128"])
        # print("citrate hydroxymutase: ", res.ssFluxesDistrib["r_4262"])
        #
        # print("malic enzyme: ", res.ssFluxesDistrib["r_0718"])
