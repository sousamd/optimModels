from geckopy import GeckoModel
from optimModels.simulation.simul_problems import GeckoSimulationProblem
from optimModels.simulation.override_simul_problem import OverrideStoicSimulProblem
from optimModels.optimization.evaluation_functions import build_evaluation_function
from optimModels.optimization.decoders import DecoderProtUnderOverExpression
import pandas
import random
from cobra.io import read_sbml_model

# LEVELS_SARA = [1e-3, 1e-2, 1e-1, 0.5, 1 , 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
LEVELS = [0, 2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5]

const_ecoli = {}


def simulate_prot():
    model = GeckoModel("single-pool")
    model.solver = 'gurobi'

    with model:
        #       for p in ["P53685","Q01574"]:
        for p in ['P33421']:
            r = model.reactions.get_by_id("draw_prot_" + p)
            r.lower_bound = 0
            r.upper_bound = 0
        res = model.optimize()
        print(" --> growth " + str(res.objective_value))
        print(" --> r_2111 " + str(res.fluxes["r_2111"]))
        print(" --> r_2056 " + str(res.fluxes["r_2056"]))
        print(" --> r_1714 " + str(res.fluxes["r_1714_REV"]))

    print(" ------------ ")


def simulate_wt(solver_name):
    model = GeckoModel('single-pool')
    model.solver = solver_name

    res = model.optimize()
    print(res)
    # for p in model.proteins:
    p = "P38066"
    with model:
        r = model.reactions.get_by_id("draw_prot_" + p)

        r.lower_bound = 0
        r.upper_bound = 0.000001
        res = model.optimize()

        # r.knock_out()
        # res = model.optimize()
        print(p + " wt simulation1 " + str(res.objective_value))

    print(str(r.lower_bound) + " --> " + str(r.upper_bound))


def simulate_wt_multi():
    some_measurements = pandas.Series({'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1})
    model = GeckoModel('multi-pool')
    model.limit_proteins(some_measurements)
    res = model.optimize()

    print(" wt simulation1 ", res.objective_value)
    for r in model.reactions:
        print(r.id, " --> ", res.fluxes[r.id])


def analysis_growth(resFileName):
    levels = [0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1]
    model = GeckoModel('single-pool')
    proteins = model.proteins
    df = pandas.DataFrame(index = proteins, columns = levels)

    for p in proteins:
        print(p)
        if p != "P38066":
            for level in levels:
                r = model.reactions.get_by_id("draw_prot_" + p)
                lb = r.lower_bound
                ub = r.upper_bound
                r.lower_bound = 0
                r.upper_bound = level
                res = model.optimize()
                df.loc[p][level] = res.objective_value
                r.lower_bound = lb
                r.upper_bound = ub
    df.to_csv(resFileName)


def analysis_ko(resFileName):
    model = GeckoModel('single-pool')
    df = pandas.DataFrame(index = range(100), columns = ["ko", "Biomass"])

    for i in range(100):
        proteins = random.sample(model.proteins, 10)
        dic = {p: 0 for p in proteins}
        model.limit_proteins(pandas.Series(dic))
        res = model.optimize()
        df.loc[i] = (dic, res.objective_value)
    df.to_csv(resFileName)


"""~~~~~~~~~~~~~~~~~~~~~~"""


def simulate_KO(model_sp, KO = [], prot_measure_fractions = None, prot_measure_ggdw = None, constraints = {},
                minbiom = 0.0, save = False):
    if model_sp == "Yeast":
        target = "r_2056"
        biom_reac = "r_2111"
        if prot_measure_fractions is None and prot_measure_ggdw is None:
            model = GeckoModel("single-pool")
        else:
            model = GeckoModel("multi-pool")
            if prot_measure_fractions:
                model.limit_proteins(fractions = prot_measure_fractions)
            else:
                model.limit_proteins(ggdw = prot_measure_ggdw)
    elif model_sp == "Ecoli":
        SBML_FILE = '../../../examples/models/eciML1515_batch.xml'
        biom_reac = 'BIOMASS_Ec_iML1515_core_75p37M'
        target = 'EX_succ_e'

        model_sbml = read_sbml_model(SBML_FILE)
        model = GeckoModel(model_sbml,
                           biomass_reaction_id = biom_reac,
                           protein_reaction_id = biom_reac,
                           carbohydrate_reaction_id = biom_reac)
    else:
        raise Exception("model_sp: Yeast or Ecoli")

    simulProb = GeckoSimulationProblem(model, constraints = constraints)

    dic = {}
    for prot in KO:
        dic["draw_prot_" + prot] = (0, 0)

    override = OverrideStoicSimulProblem(dic)
    res = simulProb.simulate(override)

    product = res.ssFluxesDistrib[target]
    biomass = res.ssFluxesDistrib[biom_reac]
    print('product', product)
    print('biomass', biomass)

    evalFunc = build_evaluation_function("WYIELD", biom_reac, target, alpha = 0.3, minBiomassValue = minbiom)
    fit, min_targ, max_targ = evalFunc.get_fitness(res, KO)
    print('fitness', fit)
    print('-------------------')

    if save:
        from Bio import SeqIO
        import urllib
        names = []
        for id in KO:
            handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/{}.xml".format(id))
            record = SeqIO.read(handle, "uniprot-xml")
            name = record.description
            names.append(name)

        robust = (min_targ / max_targ) * 100
        results_line = "KO\taero\t500\t{}\tSucc\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}%\t{}".format(model_sp, KO, biomass,
                                                                                             product, fit, max_targ,
                                                                                             min_targ, robust, names)
        save_to_csv("./results_confirm.csv", results_line)
    return res.ssFluxesDistrib, model


def simulate_UO(model_sp, UO = {}, prot_measure_fractions = None, prot_measure_ggdw = None, constraints = {},
                minbiom = 0.0, save = False):
    if model_sp == "Yeast":
        target = "r_2056"
        biom_reac = "r_2111"
        if prot_measure_fractions is None and prot_measure_ggdw is None:
            model = GeckoModel("single-pool")
        else:
            model = GeckoModel("multi-pool")
            if prot_measure_fractions:
                model.limit_proteins(fractions = prot_measure_fractions)
            else:
                model.limit_proteins(ggdw = prot_measure_ggdw)
    elif model_sp == "Ecoli":
        SBML_FILE = '../../../examples/models/eciML1515_batch.xml'
        biom_reac = 'BIOMASS_Ec_iML1515_core_75p37M'
        target = 'EX_succ_e'

        model_sbml = read_sbml_model(SBML_FILE)
        model = GeckoModel(model_sbml,
                           biomass_reaction_id = biom_reac,
                           protein_reaction_id = biom_reac,
                           carbohydrate_reaction_id = biom_reac)
    else:
        raise Exception("model_sp: Yeast or Ecoli")

    simulProb = GeckoSimulationProblem(model, constraints = constraints)

    simulProb.wt_concentrations = simulProb.simulate().get_protein_concentrations()
    criticalProteins = []
    proteins = [x for x in simulProb.model.proteins if
                x not in criticalProteins and x not in simulProb.objective.keys()]
    proteins.sort()

    decoder = DecoderProtUnderOverExpression(proteins, levels = LEVELS)

    candidate = set(decoder.decode_candidate_ids_to_index(identifiers = UO))
    override = decoder.get_override_simul_problem(candidate = candidate, simulProblem = simulProb)
    res = simulProb.simulate(override)

    product = res.ssFluxesDistrib[target]
    biomass = res.ssFluxesDistrib[biom_reac]
    print('product', product)
    print('biomass', biomass)

    evalFunc = build_evaluation_function("WYIELD", biom_reac, target, alpha = 0.3, minBiomassValue = minbiom)
    fit, min_targ, max_targ = evalFunc.get_fitness(res, UO)
    print('fitness', fit)
    print('-------------------')

    if save:
        from Bio import SeqIO
        import urllib
        names = []
        for id in UO.keys():
            handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/{}.xml".format(id))
            record = SeqIO.read(handle, "uniprot-xml")
            name = record.description
            names.append((name, UO[id][1]))

        robust = (min_targ / max_targ) * 100
        results_line = "UO\taero\t1500\t{}\tSucc\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}%\t{}".format(model_sp, UO, biomass,
                                                                                              product, fit, max_targ,
                                                                                              min_targ, robust, names)
        save_to_csv("./results_confirm.csv", results_line)
        # print(results_line)


def ecoli_model_KO(KO = [], constraints = {}):
    SBML_FILE = '../../../examples/models/eciML1515_batch.xml'
    biom_reac = 'BIOMASS_Ec_iML1515_core_75p37M'

    model = read_sbml_model(SBML_FILE)
    gecko = GeckoModel(model,
                       biomass_reaction_id = biom_reac,
                       protein_reaction_id = biom_reac,
                       carbohydrate_reaction_id = biom_reac)
    minimize = False
    simulProb = GeckoSimulationProblem(gecko, constraints = constraints, minimize = minimize)

    # WT
    res = simulProb.simulate()
    res.print()

    # ko
    dic = {}
    for prot in KO:
        dic["draw_prot_" + prot] = (0, 0)
    override = OverrideStoicSimulProblem(dic)

    res = simulProb.simulate(override)

    # print("Biotin reactions:")
    # dbts = res.ssFluxesDistrib["DBTSNo1"]
    # print('dethibiotin synth', dbts)
    # bts5 = res.ssFluxesDistrib["BTS5No1"]
    # print('biotin synth', bts5)
    # bsory = res.ssFluxesDistrib["BSORyNo1"]
    # print('btn sulf redux p', bsory)
    # bsorx = res.ssFluxesDistrib["BSORxNo1"]
    # print('btn sulf redux', bsorx)
    # btnt2ipp = res.ssFluxesDistrib["BTNt2ipp"]
    # print('btn p->c transp', btnt2ipp)
    # btnt2ipp_rev = res.ssFluxesDistrib["BTNt2ipp_REV"]
    # print('btn c->p transp', btnt2ipp_rev)
    # btntex = res.ssFluxesDistrib["BTNtex"]
    # print('btn e->p', btntex)
    # btntex_rev = res.ssFluxesDistrib["BTNtex_REV"]
    # print('btn p->e', btntex_rev)
    # btn_ex = res.ssFluxesDistrib["EX_btn_e"]
    # print('btn ex out', btn_ex)
    # btn_ex_rev = res.ssFluxesDistrib["EX_btn_e_REV"]
    # print('btn ex in', btn_ex_rev)

    print("------")
    target = 'EX_succ_e'
    target_value = res.ssFluxesDistrib[target]
    print('target: ', target_value)
    biomass = res.ssFluxesDistrib['BIOMASS_Ec_iML1515_core_75p37M']
    print('biomass: ', biomass)

    results_line = "KO\t{}\t{}\t{}\t{}\t{}\t{}".format(biomass, target, target_value, minimize, KO, constraints)
    # save_to_csv("./results.csv", results_line)
    print(results_line)


def ecoli_model_UO(UO = {}, constraints = {}):
    SBML_FILE = '../../../examples/models/eciML1515_batch.xml'
    biom_reac = 'BIOMASS_Ec_iML1515_core_75p37M'

    model = read_sbml_model(SBML_FILE)
    gecko = GeckoModel(model,
                       biomass_reaction_id = biom_reac,
                       protein_reaction_id = biom_reac,
                       carbohydrate_reaction_id = biom_reac
                       )
    minimize = False
    simulProb = GeckoSimulationProblem(gecko, constraints = constraints, minimize = minimize)

    simulProb.wt_concentrations = simulProb.simulate().get_protein_concentrations()
    criticalProteins = []
    proteins = [x for x in simulProb.model.proteins if
                x not in criticalProteins and x not in simulProb.objective.keys()]
    proteins.sort()
    LEVELS = [1e-3, 1e-2, 1e-1, 0.5, 0, 1, 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
    decoder = DecoderProtUnderOverExpression(proteins, levels = LEVELS)
    candidate = set(decoder.decode_candidate_ids_to_index(identifiers = UO))
    override = decoder.get_override_simul_problem(candidate = candidate, simulProblem = simulProb)

    res = simulProb.simulate(override)

    # print("Biotin reactions:")
    # dbts = res.ssFluxesDistrib["DBTSNo1"]
    # print('dethibiotin synth', dbts)
    # bts5 = res.ssFluxesDistrib["BTS5No1"]
    # print('biotin synth', bts5)
    # bsory = res.ssFluxesDistrib["BSORyNo1"]
    # print('btn sulf redux p', bsory)
    # bsorx = res.ssFluxesDistrib["BSORxNo1"]
    # print('btn sulf redux', bsorx)
    # btnt2ipp = res.ssFluxesDistrib["BTNt2ipp"]
    # print('btn p->c transp', btnt2ipp)
    # btnt2ipp_rev = res.ssFluxesDistrib["BTNt2ipp_REV"]
    # print('btn c->p transp', btnt2ipp_rev)
    # btntex = res.ssFluxesDistrib["BTNtex"]
    # print('btn e->p', btntex)
    # btntex_rev = res.ssFluxesDistrib["BTNtex_REV"]
    # print('btn p->e', btntex_rev)
    # btn_ex = res.ssFluxesDistrib["EX_btn_e"]
    # print('btn ex out', btn_ex)
    # btn_ex_rev = res.ssFluxesDistrib["EX_btn_e_REV"]
    # print('btn ex in', btn_ex_rev)

    print("------")
    target = 'EX_succ_e'
    target_value = res.ssFluxesDistrib[target]
    print('target: ', target_value)
    biomass = res.ssFluxesDistrib['BIOMASS_Ec_iML1515_core_75p37M']
    print('biomass: ', biomass)

    results_line = "UO\t{}\t{}\t{}\t{}\t{}\t{}".format(biomass, target, target_value, minimize, UO, constraints)
    save_to_csv("./results.csv", results_line)


def save_to_csv(csv_address, line):
    file = open(csv_address, "a+")
    file.write("\n" + line)
    file.close()


if __name__ == "__main__":
    # SBML_FILE = '../../../examples/models/eciML1515_batch.xml'
    # biom_reac = 'BIOMASS_Ec_iML1515_core_75p37M'
    #
    # model = read_sbml_model(SBML_FILE)
    # gecko = GeckoModel(model,
    #                    biomass_reaction_id = biom_reac,
    #                    protein_reaction_id = biom_reac,
    #                    carbohydrate_reaction_id = biom_reac)
    #
    # simulProb = GeckoSimulationProblem(gecko, constraints=const_ecoli)
    #
    #
    # essenciais = simulProb.find_essential_proteins()
    #
    # print(essenciais)
    #
    proteins1 = ['P0A725', 'P0A9L8', 'P0AC75', 'P0ACB2', 'P21156', 'P38038', 'P17846',
                 'P0A786', 'P0A9D8', 'P00861', 'P07004', 'P00934', 'P0AAI9', 'P08244', 'P0A8K1',
                 'P0AEA8', 'P07623', 'P61714', 'P30126', 'P27300', 'P23908', 'P0AB89', 'P23845',
                 'P0AGK1', 'P0A6C5', 'P0ACB4', 'P60757', 'P0AF98', 'P0AEZ1', 'P10371', 'P0ADC6',
                 'P0A6A6', 'P0ABD5', 'P0A6L2', 'P45578', 'P06989', 'P60546', 'P0ACC7', 'P05194',
                 'P05791', 'P04036', 'P0AC13', 'P11458', 'P00547', 'P0ABD8', 'P22634', 'P0A9D4',
                 'P0A6C8', 'P06988', 'P0A6E4', 'P15770', 'P08200', 'P0AG40', 'P62615', 'P23830',
                 'P0A7D4', 'P0AD68', 'P0A7E3', 'P0A6I9', 'P0AEK2', 'P22255', 'P0A7I7', 'P62617',
                 'P26281', 'P0ADC1', 'P23893', 'P08373', 'P11447', 'P19624', 'P04079', 'P0ADV9',
                 'P0A6I6', 'P29680', 'P0A9V1', 'P31057', 'P0A6I3', 'P09029', 'P0A9Q5', 'P00909',
                 'P05020', 'P0A6A8', 'P11880', 'P0A9Q9', 'P04951', 'P0A6J1', 'P17854', 'P0A790',
                 'P24182', 'P0AED7', 'P17169', 'P0ABZ4', 'P0AD65', 'P17952', 'P09151', 'P43341',
                 'P0AB80', 'P08192', 'P25539', 'P22939', 'P07639', 'P0ABG4', 'P0A749', 'P0AF12',
                 'P05793', 'P31663', 'P60752', 'P0AGG0', 'P11446', 'P00904', 'P0A7B3', 'P77488',
                 'P00903', 'P17443', 'Q46893', 'P60664', 'P0A9J8', 'P0A7B5', 'P0A884', 'P62620',
                 'P07023', 'P06983', 'P0A6T5', 'P0A794', 'P60595', 'P0AF16', 'P0A7J0', 'P05041',
                 'P06986', 'P12008', 'P0A6K1', 'P00935', 'P0A879', 'P31120', 'P0A6X1', 'P0A715',
                 'P0A6D3', 'P60472', 'P10441', 'P30011', 'P0A752', 'P0ABG1', 'P0ABQ0', 'P0A7E5',
                 'P22188', 'P0A6W3', 'P31554', 'P62623', 'P04805', 'P21645', 'P10902', 'P14900',
                 'P0AFU8', 'P26647', 'P0ADV1', 'P0A722', 'P0A877', 'P0A817', 'P45568', 'P06987',
                 'P18843', 'P0A720', 'P00895', 'P0AD57']
    proteins2 = ['P00934', 'P0ABZ4', 'P17846', 'P0ACB2', 'P0A8K1', 'P0AEK2', 'P0AGK1', 'P0A7E5',
                 'P25539', 'P0ADV1', 'P0A817', 'P0A7J0', 'P0A9Q9', 'P0ADC6', 'P11880', 'P17443', 'P0A9Q5',
                 'P0A749', 'P21156', 'P0AF12', 'P0A6I6', 'Q46893', 'P30011', 'P04079', 'P0A794', 'P08373',
                 'P60595', 'P04951', 'P0AB89', 'P0A6E4', 'P21645', 'P0A6L2', 'P0A752', 'P60752', 'P0AEZ1',
                 'P0A720', 'P31663', 'P0AC75', 'P0A6C5', 'P04036', 'P05020', 'P23908', 'P0ACC7', 'P06988',
                 'P06983', 'P45568', 'P09151', 'P0A6X1', 'P0AAI9', 'P14900', 'P0A7I7', 'P0A6T5', 'P0AF16',
                 'P0AFU8', 'P06986', 'P24182', 'P0ABH7', 'P00547', 'P0AD68', 'P0A6J1', 'P17952', 'P0A6A6',
                 'P0ABQ0', 'P22634', 'P0A6I3', 'P00935', 'P22188', 'P23830', 'P0A6K1', 'P38038', 'P60664',
                 'P11458', 'P0A6A8', 'P0ABG1', 'P00904', 'P07023', 'P0A790', 'P07639', 'P11447', 'P0AC13',
                 'P0A722', 'P0AG40', 'P08244', 'P0A879', 'P0A9J8', 'P0ADC1', 'P0A6C8', 'P23893', 'P26647',
                 'P18843', 'P17169', 'P0ADV9', 'P07623', 'P00903', 'P77488', 'P30126', 'P00909', 'P17854',
                 'P61714', 'P62620', 'P0A715', 'P0AD65', 'P09126', 'P11446', 'P08192', 'P31120', 'P26281',
                 'P29680', 'P12008', 'P60757', 'P0AED7', 'P60472', 'P0A6W3', 'P05793', 'P0ABG4', 'P06987',
                 'P0A786', 'P0A7B5', 'P0ABD5', 'P43341', 'P31554', 'P0A877', 'P0A6D3', 'P05791', 'P05041',
                 'P62623', 'P45578', 'P22939', 'P19624', 'P04805', 'P23845', 'P60546', 'P05194', 'P0A725',
                 'P0A9V1', 'P0AEA8', 'P10441', 'P0ABD8', 'P62615', 'P0AD57', 'P10902', 'P00861', 'P0AF98',
                 'P31057', 'P0AB80', 'P0A7E3', 'P0AGG0', 'P07004', 'P0A6I9', 'P0AC16', 'P0A884', 'P0A7B3',
                 'P00895', 'P62617', 'P0ACB4', 'P06989', 'P27300', 'P10371', 'P15770', 'P0A9L8', 'P08200', 'P0A9D8',
                 'P0A9D4']
    #
    # print([x for x in proteins1 if x in proteins2])

    # Yeast
    const_aero = {'r_1714_REV': (0, 10)}  # glucose
    const_anaerobic = {
        'r_1714_REV': (0, 10),  # glucose
        'r_1992_REV': (0, 0),  # oxygen
        'r_1994_REV': (0, 1000), 'r_2189_REV': (0, 1000), 'r_1757_REV': (0, 1000)}  # lipids

    min_biom = 0.03135  # aero
    # min_biom = 0.01875  #anaero

    # # Ecoli
    # const_aero = {'EX_glc__D_e_REV': (0, 10)}
    # const_anaero = {'EX_glc__D_e_REV': (0, 10), 'EX_o2_e_REV': (0, 0)}
    #
    # min_biom = 0.03112  # aero
    # # min_biom = 0.00518  # anaero

    KO = ['P42941', 'P54115', 'P37299', 'P37303', 'P46969', 'Q12680', 'P14065']
    UO = {'P60752': (0, 2), 'P00864': (0, 32), 'P0ABB4': (0, 32), 'P0AD65': (0, 0.25), 'P46022': (0, 32),
          'P16701': (0, 0.5), 'P0ABG1': (0, 8), 'P06983': (0, 32), 'P05042': (0, 0.25), 'P0ABP8': (0, 0.0625)}

    # saveflag = True
    saveflag = False

    fluxes, model = simulate_KO(model_sp = "Yeast", KO = KO, constraints = const_aero, minbiom = min_biom,
                                save = saveflag)
    # simulate_UO(model_sp = "Ecoli", UO = UO, constraints = const_aero, minbiom = min_biom, save = saveflag)
