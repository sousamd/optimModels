from optimModels.model.kineticModel import load_kinetic_model
from framed.cobra.simulation import FBA, pFBA
from framed.io.sbml import load_cbmodel
from optimModels.simulation.simul_problems import StoicSimulationProblem, KineticSimulationProblem
from optimModels.simulation.override_simul_problem import  OverrideKineticSimulProblem
from optimModels.utils.utils import fix_exchange_reactions_model
from cobra.io import read_sbml_model
from optimModels.simulation.override_simul_problem import OverrideStoicSimulProblem
from optimModels.optimization.evaluation_functions import build_evaluation_function
from framed import essential_reactions

##################################################################################
# File used to validate some of the results obtained in the optimization process #
##################################################################################
from optimModels.unittests.cbm_optim import get_nontargets

# Simulation:
# model: comunity model
# objective: Maximize succinate production with a biomass >0.55
def medium_test():
    SBML_FILE = "../../../examples/models/EC_SC_model.xml"
    model = load_cbmodel(SBML_FILE)
    model2 = fix_exchange_reactions_model(model)

    for r in model2.reactions:
        print (r)

    for r_id, rxn in model2.reactions.items():
        if r_id.startswith('R_EX_'):  # ou tambem podes fazer if rxn.is_exchange:
            rxn.lb = 0
            rxn.ub = 1000 \
                #if rxn.ub is None else rxn.ub


    #sol 1
    constraints = {
    'R_EX_mn2_e__mod1':(-1000,0), 'R_EX_cobalt2_e__mod1':(-1000,0), 'R_EX_cl_e__mod1':(-1000,0), 'R_EX_xyl_D_medium_':(-12,0),
    'R_EX_mobd_e__mod1':(-1000,0), 'R_EX_mg2_e__mod1':(-1000,0),
    'R_EX_o2_medium_':(-2,0), 'R_EX_k_medium_':(-1000,0),
    'R_EX_nh4_medium_':(-1000,0), 'R_EX_fe3_e__mod1':(-1000,0),
    'R_EX_zn2_e__mod1':(-1000,0), 'R_EX_ca2_e__mod1':(-1000,0), 'R_EX_so4_medium_':(-1000,0), 'R_EX_pi_medium_':(-1000,0),
    'R_EX_cu2_e__mod1':(-1000,0)}

    #sol2
    constraints = {
    'R_EX_mn2_e__mod1':(-1000,0), 'R_EX_cobalt2_e__mod1':(-1000,0), 'R_EX_cl_e__mod1':(-1000,0), 'R_EX_mobd_e__mod1':(-1000,0), 'R_EX_cu_e__mod1':(-1000,0),
    'R_EX_mg2_e__mod1':(-1000,0), 'R_EX_glc_D_medium_':(-12,0), 'R_EX_o2_medium_':(-2,0), 'R_EX_k_medium_':(-1000,0), 'R_EX_nh4_medium_':(-1000,0),
    'R_EX_fe3_e__mod1':(-1000,0), 'R_EX_zn2_e__mod1':(-1000,0), 'R_EX_ca2_e__mod1':(-1000,0), 'R_EX_so4_medium_':(-1000,0), 'R_EX_pi_medium_':(-1000,0)}

    res  = pFBA(model2, objective={"R_BM_total_Synth": 1}, constraints=constraints)
    #res  = FBA(model2, objective={"R_BM_total_Synth": 1})

    for r, f in res.values.items():
        print (r + " --> " +str(f))
    print(res.values["R_EX_succ_medium_"])
    print(res.values["R_BM_total_Synth"])
    print(res.values["R_EX_glc_D_medium_"])

# Simulation:
# model: Ec_iAF1260
# objective: Maximize succinate production with a biomass >0.55
def cbm_simualtion_iAF():
    SBML_FILE = "../../../examples/models/Ec_iAF1260.xml"
    #model = load_cbmodel(SBML_FILE, flavor="cobra")
    model = read_sbml_model(SBML_FILE)
    #constraints = {'R_Ec_biomass_iAF1260_core_59p81M': (0.55, 9999)}
    '''
    res  = FBA(model, objective={"R_EX_succ_e": 1}, constraints=constraints)

    for r, f in res.values.items():
        if f !=0:
            print (r + " --> " +str(f))

    print(res.values["R_EX_succ_e"])
    print(res.values["R_Ec_biomass_iAF1260_core_59p81M"])
    '''
    simulProblem = StoicSimulationProblem(model, objective={"Ec_biomass_iAF1260_core_59p81M": 1}, constraints=None, withCobraPy=True)
    from optimModels.simulation.override_simul_problem import OverrideStoicSimulProblem
    listas_new = [['GPDDA3'],
                  ['3OAR140', 'KG6PDC', 'NOtpp', 'DURIt2pp'],
                  ['G3PItex', 'ACGS', 'PGI', 'G3PCabcpp', 'D_LACt2pp'],
                  ['2AGPGAT161', 'EX_imp_e_', 'GPDDA5pp', 'G3PAT141'],
                  ['FACOAE120'],
                  ['THMDt2pp', 'PPK2r', 'PLIPA2G141pp', 'EX_xmp_e_'],
                  ['ORNabcpp', 'SULFACtex', 'PYAM5PO', '2AGPA140tipp', 'ACBIPGT'],
                  ['EX_alaala_e_', 'MSO3tex', '2AGPE180tipp', 'EX_pser_L_e_'],
                  ['OCTAtex'],
                  ['HCINNMt2rpp'],
                  ['HDCOAI', 'GALTptspp'],
                  ['EX_gly_e_', 'UPPDC1', 'EX_gmp_e_', 'EX_feoxam_un_e_'],
                  ['ADK1', 'CRNDCAL2', 'Ktex', 'GAL1PPpp'],
                  ['ADK1', 'CRNDCAL2', 'Ktex', 'GAL1PPpp']]

    dics = []
    for lista in listas_new:
        dic = {}
        for reac in lista:
            dic[reac] = (0, 0)
        dics.append(dic)
    import time
    start = time.time()
    for dic in dics:
        simulProblem = StoicSimulationProblem(model, objective={"Ec_biomass_iAF1260_core_59p81M": 1}, constraints=None, withCobraPy=True)
        override = OverrideStoicSimulProblem(dic)
        res = simulProblem.simulate(override)
        print(res.ssFluxesDistrib['EX_succ_e_'])
        #model = model.copy()
        model = read_sbml_model(SBML_FILE)
    print("--- %s seconds 1core ---" % (time.time() - start), 'end')
    #print(model.reactions.r_0630.bounds)


# Simulation:
# model: iMM904
# objective: Ethanol production with medium definition (other uptake reactions are defined in the model) and limit of biomass production
def cbm_simualtion():
    SBML_FILE = "../../../examples/models/iMM904.xml"
    model = load_cbmodel(SBML_FILE, flavor="fbc2")

    constraints = {'R_EX_so4_e': (-10000,0), 'R_EX_o2_e': (-50,0), 'R_EX_gam6p_e': (-10,0), 'R_EX_melib_e': (-5,0)}
    constraints.update({"R_BIOMASS_SC5_notrace":(0.21,9999)})
    res  = FBA(model, objective={"R_EX_etoh_e": 1}, constraints=constraints)
    print(res.values["R_EX_etoh_e"])
    print(res.values["R_BIOMASS_SC5_notrace"])
    print(res.values["R_EX_so4_e"])
    print(res.values["R_EX_gam6p_e"])
    print(res.values["R_EX_o2_e"])
    print(res.values["R_EX_melib_e"])


def simul_iJO1366SL(KO, constraints=None, withCobraPy=False):
    SBML_FILE = '../../../examples/models/iJO1366SL.xml'
    model = load_cbmodel(SBML_FILE, flavor='cobra')
    model.set_reaction_objective('R_Ec_biomass_iJO1366_core_53p95M', coeff=1)

    criticalReacs = get_nontargets('../../../examples/models/nontargets_iJO1366SL_UO.txt')

    simulProblem = StoicSimulationProblem(model, constraints=constraints)

    LEVELS = [1e-3, 1e-2, 1e-1, 0.5, 1, 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
    reactions = [x for x in simulProblem.get_internal_reactions() if x not in criticalReacs and x not in simulProblem.objective.keys()]
    reactions.sort()


    from optimModels.optimization.decoders import DecoderReacUnderOverExpression
    decoder = DecoderReacUnderOverExpression(reactions, levels=LEVELS)


    candidate = decoder.decode_candidate_ids_to_index(identifiers=KO)
    print(candidate)

    override = decoder.get_override_simul_problem(candidate=candidate, simulProblem=simulProblem)

    #res_wt = simulProblem.simulate()
    #print(res_wt.ssFluxesDistrib['R_GARFT'])

    res = simulProblem.simulate(override)
    res.print()
    print(res.ssFluxesDistrib['R_EX_succ_e'])
    print(res.ssFluxesDistrib['R_Ec_biomass_iJO1366_core_53p95M'])
    print((res.ssFluxesDistrib['R_EX_succ_e']*res.ssFluxesDistrib['R_Ec_biomass_iJO1366_core_53p95M'])/(-res.ssFluxesDistrib['R_EX_glc_e']))
    #print(res_wt.ssFluxesDistrib['R_GARFT'])


    #print(reactions.index('R_ACKr'))


    # from framed import essential_reactions
    # essential_reacs = essential_reactions(model)
    # print(essential_reacs)

    # for rid, robj in model.reactions.items():
    #     if robj.lb is None: robj.lb = -100000.0
    #     if robj.ub is None: robj.ub = 100000.0
    #
    # #print(model.reactions.R_TALA.lb,model.reactions.R_TALA.ub)
    #
    #simulProblem_wt = StoicSimulationProblem(model, method='FBA')
    #
    # dic = {}
    # for reac in KO:
    #     dic[reac] = (0,0)
    #
    #wt = simulProblem_wt.simulate()
    #fwt = wt.ssFluxesDistrib['R_PGI']
    #wt.print()
    #
    # simulProblem = StoicSimulationProblem(model)
    # dic = {'R_PGI':(fwt*2, fwt*2)}
    #
    # override = OverrideStoicSimulProblem(dic)
    # res = simulProblem.simulate(override)
    # res.print()
    # print(res.ssFluxesDistrib['R_EX_succ_e'])
    # print(res.ssFluxesDistrib['R_Ec_biomass_iJO1366_core_53p95M'])
    # print((res.ssFluxesDistrib['R_EX_succ_e']*res.ssFluxesDistrib['R_Ec_biomass_iJO1366_core_53p95M'])/(-res.ssFluxesDistrib['R_EX_glc_e']))

def simul_iMM904SL(KO, constraints=None, withCobraPy=False):
    SBML_FILE = '../../../examples/models/iMM904SL_v6.xml'
    model = load_cbmodel(SBML_FILE, flavor='cobra')
    model.set_reaction_objective('R_biomass_SC5_notrace', coeff=1)

    model_cobra = read_sbml_model(SBML_FILE)

    criticalReacs = get_nontargets('../../../examples/models/nontargets_iMM904SL_v6_UO.txt')

    simulProblem = StoicSimulationProblem(model_cobra, constraints=constraints, withCobraPy=withCobraPy)



    #LEVELS = [1e-3, 1e-2, 1e-1, 0.5, 0, 1, 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
    reactions = [x for x in simulProblem.get_internal_reactions() if
                 x not in criticalReacs and x not in simulProblem.objective.keys()]
    reactions.sort()

    # from optimModels.optimization.decoders import DecoderReacUnderOverExpression, DecoderReacKnockouts
    # decoder = DecoderReacUnderOverExpression(reactions, levels=LEVELS)
    # #decoder = DecoderReacKnockouts(reactions)
    #
    # candidate = decoder.decode_candidate_ids_to_index(identifiers=KO)
    # override = decoder.get_override_simul_problem(candidate=candidate, simulProblem=simulProblem)


    res = simulProblem.simulate()

    print(res.ssFluxesDistrib['R_EX_succ_e'])
    print(res.ssFluxesDistrib['R_biomass_SC5_notrace'])

    # print((res.ssFluxesDistrib['R_EX_succ_e'] * res.ssFluxesDistrib['R_biomass_SC5_notrace']) / (
    #       -res.ssFluxesDistrib['R_EX_glc_e']))



def simplify_sol_iJO1366SL(candidate, KO=True, func="BPCY", constraints=None):
    SBML_FILE = '../../../examples/models/iJO1366SL.xml'
    model = load_cbmodel(SBML_FILE, flavor='cobra')
    model.set_reaction_objective('R_Ec_biomass_iJO1366_core_53p95M', coeff=1)
    simulProb = StoicSimulationProblem(model, constraints=constraints)
    if func == "BPCY":
        objFunction = build_evaluation_function("BPCY", "R_Ec_biomass_iJO1366_core_53p95M", "R_EX_succ_e", "R_EX_glc_e")
    else:
        objFunction = build_evaluation_function("targetFlux", ["R_EX_succ_e"])

    dic = {}
    if KO:
        for prot in candidate:
            dic[prot] = (0, 0)
    else:
        for prot in candidate.keys():
            dic[prot] = candidate[prot]

    override_all = OverrideStoicSimulProblem(dic)
    res_all = simulProb.simulate(override_all)
    print(res_all.ssFluxesDistrib['R_EX_succ_e'])
    print(res_all.ssFluxesDistrib['R_Ec_biomass_iJO1366_core_53p95M'])
    fitness = (res_all.ssFluxesDistrib['R_EX_succ_e']* res_all.ssFluxesDistrib['R_Ec_biomass_iJO1366_core_53p95M'])/-res_all.ssFluxesDistrib['R_EX_glc_e']
    print(fitness)

    new_dic = dic.copy()
    for sol in dic:
        del new_dic[sol]
        override = OverrideStoicSimulProblem(new_dic)
        res = simulProb.simulate(override)
        newFitness = objFunction.get_fitness(res, override)
        #print(newFitness)
        if round(fitness, 5) != round(newFitness, 5):
            new_dic[sol] = dic[sol]
    if KO:
        prots = []
        for key in new_dic:
            if key not in constraints.keys(): prots.append(key)
        return prots
    else:
        x = {}
        for key in new_dic:
            if key not in constraints.keys(): x[key] = new_dic[key]
        return x


def simplify_sol_iMM904SL(candidate, KO=True, func="BPCY", constraints=None):
    SBML_FILE = '../../../examples/models/iMM904SL_v6.xml'
    model = load_cbmodel(SBML_FILE, flavor='cobra')
    model.set_reaction_objective('R_biomass_SC5_notrace', coeff=1)
    simulProb = StoicSimulationProblem(model, constraints=constraints)
    if func == "BPCY":
        objFunction = build_evaluation_function("BPCY", "R_biomass_SC5_notrace", "R_EX_succ_e", "R_EX_glc_e")
    else:
        objFunction = build_evaluation_function("targetFlux", ["r_2056"])
    dic = {}

    if KO:
        for prot in candidate:
            dic[prot] = (0, 0)
    else:
        for prot in candidate.keys():
            dic[prot] = candidate[prot]


    override_all = OverrideStoicSimulProblem(dic)
    res_all = simulProb.simulate(override_all)
    print(res_all.ssFluxesDistrib['R_EX_succ_e'])
    print(res_all.ssFluxesDistrib['R_biomass_SC5_notrace'])
    fitness = (res_all.ssFluxesDistrib['R_EX_succ_e']* res_all.ssFluxesDistrib['R_biomass_SC5_notrace'])/-res_all.ssFluxesDistrib['R_EX_glc_e']
    print(fitness)

    new_dic = dic.copy()
    for sol in dic:
        del new_dic[sol]
        override = OverrideStoicSimulProblem(new_dic)
        res = simulProb.simulate(override)
        newFitness = objFunction.get_fitness(res, override)
        if round(fitness, 8) != round(newFitness, 8):
            new_dic[sol] = dic[sol]
    if KO:
        prots = []
        for key in new_dic:
            if key not in constraints.keys(): prots.append(key)
        return prots
    else:
        x = {}
        for key in new_dic:
            if key not in constraints.keys(): x[key] = new_dic[key]
        return x




if __name__ == '__main__':
    #ECOLI
    #const_ecoli = {'R_EX_glc_e': (-10.0, 100000.0), 'R_EX_o2_e': (-9.66, 100000.0)}
    #simul_iJO1366SL(KO=UO, constraints=const_ecoli)
    #print(simplify_sol_iJO1366SL(candidate=UO, constraints=const_ecoli, KO=False))

    #YEAST
    const_yeast = {'R_EX_glc_e': (-10.0, 999999.0), 'R_EX_o2_e': (0,0)} #'R_EX_o2_e': (-12.25, 100000.0)}
    #
    # lista = {'R_ALCD2x':(0,0),
    #          'R_ALCD19y':(0,0),'R_G3PD1ir':(0,0)}

    simul_iMM904SL(KO=[], constraints=const_yeast, withCobraPy=True)

    #res = simplify_sol_iMM904SL(candidate=lista, constraints= const_yeast, KO=False)
    #print(res)
