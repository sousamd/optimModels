from geckopy import GeckoModel
from optimModels.simulation.simul_problems import GeckoSimulationProblem
from optimModels.simulation.override_simul_problem import OverrideStoicSimulProblem
from optimModels.optimization.evaluation_functions import build_evaluation_function
from optimModels.optimization.decoders import DecoderProtUnderOverExpression, DecoderProtKnockouts
import pandas
import random
from cobra.io import read_sbml_model

# LEVELS_SARA = [1e-3, 1e-2, 1e-1, 0.5, 1 , 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
LEVELS = [0, 2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5, 0.1, 10]


def simulate_prot():
    model = GeckoModel("single-pool")
    model.solver = 'cplex'

    with model:
        # for p in ["P53685","Q01574"]:
        for p in ['P33421']:
            r = model.reactions.get_by_id("draw_prot_" + p)
            r.lower_bound = 0
            r.upper_bound = 0
        res = model.optimize()
        print(" --> growth "+ str(res.objective_value))
        print (" --> r_2111 "+ str(res.fluxes["r_2111"]))
        print(" --> r_2056 " + str(res.fluxes["r_2056"]))
        print (" --> r_1714 "+ str(res.fluxes["r_1714_REV"]))

    print(" ------------ ")

def simulate_wt ():
    model = GeckoModel('single-pool')
    res = model.optimize()
    print (res)
    #for p in model.proteins:
    p= "P38066"
    with model:
        r = model.reactions.get_by_id("draw_prot_" + p)

        lb = r.lower_bound
        ub = r.upper_bound
        r.lower_bound = 0
        r.upper_bound = 0.000001
        res = model.optimize()

        #r.knock_out()
        #res = model.optimize()
        print( p + " wt simulation1 " + str(res.objective_value))

    print (str(r.lower_bound) +" --> " +  str(r.upper_bound))

def simulate_wt_multi ():
    model = GeckoModel('multi-pool')
    import pandas
    some_measurements = pandas.Series({'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1})
    model = GeckoModel('multi-pool')
    model.limit_proteins(some_measurements)
    res = model.optimize()

    print(" wt simulation1 " , res.objective_value)
    for r in model.reactions:
        print (r.id," --> " ,res.fluxes[r.id])

def analysis_growth (resFileName):
    levels = [0,1e-10,1e-8,1e-6,1e-4,1e-2,0.1]
    model = GeckoModel('single-pool')
    proteins = model.proteins
    df = pandas.DataFrame(index=proteins, columns=levels)

    for p in proteins:
        print(p)
        if p !="P38066":
            for level in levels:
                r = model.reactions.get_by_id("draw_prot_" + p)
                lb = r.lower_bound
                ub = r.upper_bound
                r.lower_bound = 0
                r.upper_bound = level
                res = model.optimize()
                df.loc[p][level]= res.objective_value
                r.lower_bound = lb
                r.upper_bound = ub
    df.to_csv(resFileName)

def analysis_ko (resFileName):
    levels = [0,1e-10,1e-8,1e-6,1e-4,1e-2,0.1]

    model =  GeckoModel('single-pool')
    df = pandas.DataFrame(index = range(100), columns=["ko", "Biomass"])

    for i in range(100):
        proteins = random.sample(model.proteins,10)
        dic = {p:0 for p in proteins}
        model.limit_proteins(pandas.Series(dic))
        res = model.optimize()
        df.loc[i]= (dic, res.objective_value)
    df.to_csv(resFileName)

def simulate_KO(target, KOs=[[]], prot_measure_fractions=None, prot_measure_ggdw=None, constraints={}, minBiomassValue=None):
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)

    simulProb = GeckoSimulationProblem(model, constraints=constraints)

    fits = []
    bios = []
    succs_bio = []
    maxs = []
    mins = []
    solutions = []
    for candidate in KOs:
        dic = {}
        for prot in candidate:
            if "draw_prot_" + prot in model.reactions:
                dic["draw_prot_" + prot] = (0, 0)
            else:
                dic["prot_" + prot + '_exchange'] = (0, 0)

        override = OverrideStoicSimulProblem(dic)
        res = simulProb.simulate(override)

        succ_bio = res.ssFluxesDistrib[target]
        bio = res.ssFluxesDistrib['r_2111']

        print('product', succ_bio)
        print('biomass', bio)

        bios.append(bio)
        succs_bio.append(succ_bio)


        evalFunc = build_evaluation_function("WYIELD", "r_2111", target, alpha=0.3, minBiomassValue=minBiomassValue)
        fit = evalFunc.get_fitness(res, candidate)
        print('fitness', fit)
        print('-------------------')

        fits.append(fit[0])
        maxs.append(fit[1])
        mins.append(fit[2])

        uniprot = ''
        for prot in candidate:
            prot_name = get_uniprot_name(prot)
            uniprot += prot + ':' + prot_name + ', '
        solutions.append(uniprot)



    df = pandas.DataFrame({'biomass':bios, 'fitness':fits, 'succ_bio':succs_bio,  'max_succ':maxs, 'min_succ':mins,
                           'candidate':solutions}, columns=['fitness','biomass', 'succ_bio', 'max_succ', 'min_succ', 'candidate'])
    df.to_excel('C:/Users/BiSBII/Results/multi_pool/aero_multi/simulation.xlsx')

def simulate_UO(target, wt_concentrations=[], wt_fluxes=[], UOs=[], prot_measure_fractions=None, prot_measure_ggdw=None, constraints={}, minBiomassValue=None, objective={'r_2111':1}, minimize=False):

    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)


    simulProb = GeckoSimulationProblem(model, constraints=constraints, objective=objective, minimize=minimize)


    simulProb.wt_concentrations = simulProb.simulate().get_protein_concentrations()
    simulProb.wt_fluxes = wt_fluxes
    simulProb._protein_rev_reactions = simulProb.get_protein_rev_reactions()


    criticalProteins = []
    proteins = [x for x in simulProb.model.proteins if x not in criticalProteins and x not in simulProb.objective.keys()]
    proteins.sort()

    decoder = DecoderProtUnderOverExpression(proteins, levels=LEVELS)


    fits = []
    bios = []
    succs_bio = []
    maxs = []
    mins = []
    solutions = []

    for UO in UOs:
        candidate = set(decoder.decode_candidate_ids_to_index(identifiers=UO))
        override = decoder.get_override_simul_problem(candidate=candidate, simulProblem=simulProb)
        res = simulProb.simulate(override)
        succ_bio = res.ssFluxesDistrib[target]
        bio = res.ssFluxesDistrib['r_2111']



        print('product', succ_bio)
        print('biomass', bio)
        #print('glucose', res.ssFluxesDistrib['r_1714_REV'])
        print('DIR')
        print(res.ssFluxesDistrib['r_0659No1'])

        print('REV')
        print(res.ssFluxesDistrib['r_0659_REVNo1'])


        # bios.append(bio)
        # succs_bio.append(succ_bio)



        # evalFunc = build_evaluation_function("WYIELD", "r_2111", "r_2056", alpha=0.3, minBiomassValue=minBiomassValue)
        # fit = evalFunc.get_fitness(res, UO)
        # print('fitness', fit)
        # print('-------------------')
    #
    #     fits.append(fit[0])
    #     maxs.append(fit[1])
    #     mins.append(fit[2])
    #
    #     uniprot = ''
    #     for prot in UO:
    #         prot_name = get_uniprot_name(prot)
    #         uniprot += prot + '(' + prot_name + '):' + str(UO[prot])+ ', '
    #     solutions.append(uniprot)
    #
    # df = pandas.DataFrame(
    #          {'biomass': bios, 'fitness': fits, 'succ_bio': succs_bio, 'max_succ': maxs, 'min_succ': mins,
    #           'candidate': solutions}, columns=['fitness', 'biomass', 'succ_bio', 'max_succ', 'min_succ', 'candidate'])
    # df.to_excel('C:/Users/BiSBII/Results/multi_pool/aero_multi/simulation.xlsx')

def verify_sol(solution, target, constraints={}, prot_measure_fractions =None, prot_measure_ggdw= None, KO=True, func="BPCY", minBiomassValue=None):
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)

    simulProb = GeckoSimulationProblem(model, constraints=constraints)
    if func == "BPCY":
        objFunction = build_evaluation_function("BPCY", "r_2111", target, "r_1714_REV")
    elif func == 'tf':
        objFunction = build_evaluation_function("targetFlux", [target])
    elif func == 'wyield':
        objFunction = build_evaluation_function("WYIELD", "r_2111", target, alpha=0.3, minBiomassValue=minBiomassValue)

    dic = {}
    if KO:
        for prot in solution:
            dic["draw_prot_" + prot] = (0, 0)

        override_all = OverrideStoicSimulProblem(dic)
        res_all = simulProb.simulate(override_all)
        fitness = objFunction.get_fitness(res_all, override_all)[0]
        print(res_all.ssFluxesDistrib[target], "product")
        print(res_all.ssFluxesDistrib['r_2111'], "biomass")
        print(fitness)

        new_dic = dic.copy()
        for sol in dic:
            if sol not in constraints:
                del new_dic[sol]
                override = OverrideStoicSimulProblem(new_dic)
                res = simulProb.simulate(override)
                newFitness = objFunction.get_fitness(res, override)[0]
                print('fit', newFitness, 'prot', sol)
                print('____________')
                if round(fitness, 5) != round(newFitness, 5):
                    new_dic[sol] = dic[sol]


        prots = []
        for key in new_dic:
            if key not in constraints.keys(): prots.append(key[10:])
        return prots

    else:
        simulProb.wt_concentrations = simulProb.simulate().get_protein_concentrations()
        criticalProteins = []
        proteins = [x for x in simulProb.model.proteins if
                    x not in criticalProteins and x not in simulProb.objective.keys()]
        proteins.sort()
        #LEVELS = [1e-3, 1e-2, 1e-1, 0.5, 0, 1, 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
        decoder = DecoderProtUnderOverExpression(proteins, levels=LEVELS)
        candidate = set(decoder.decode_candidate_ids_to_index(identifiers=solution))
        override_all = decoder.get_override_simul_problem(candidate=candidate, simulProblem=simulProb)
        res_all = simulProb.simulate(override_all)
        fitness = objFunction.get_fitness(res_all, override_all)[0]
        print(res_all.ssFluxesDistrib[target], "product")
        print(res_all.ssFluxesDistrib['r_2111'], "biomass")
        print(fitness)


        new_cand = list(candidate)
        for sol in list(candidate):
            new_cand.remove(sol)
            override = decoder.get_override_simul_problem(candidate=new_cand, simulProblem=simulProb)
            res = simulProb.simulate(override)
            newFitness = objFunction.get_fitness(res, override)[0]
            print('fit', newFitness, 'prot', proteins[sol[0]])
            print('____________')
            if round(fitness, 5) != round(newFitness, 5):
                new_cand.append(sol)


        final_override = decoder.get_override_simul_problem(new_cand, simulProb)
        return final_override.get_modifications()


def ecoli_model_KO(KO=[], constraints={}):
    SBML_FILE = '../../../examples/models/eciML1515_batch.xml'

    model = read_sbml_model(SBML_FILE)
    gecko = GeckoModel(model, biomass_reaction_id='BIOMASS_Ec_iML1515_core_75p37M')
    simulProb = GeckoSimulationProblem(gecko, constraints=constraints)

    essenciais = simulProb.find_essential_proteins()
    print(essenciais)
    '''
    dic = {}
    for prot in KO:
        dic["draw_prot_" + prot] = (0, 0)
    override = OverrideStoicSimulProblem(dic)
    res = simulProb.simulate(override)
    print('product', res.ssFluxesDistrib['EX_succ_e'])
    print('biomass', res.ssFluxesDistrib['BIOMASS_Ec_iML1515_core_75p37M'])
    '''

def ecoli_model_UO(UO={}, constraints={}):
    SBML_FILE = '../../../examples/models/eciML1515_batch.xml'
    model = read_sbml_model(SBML_FILE)
    gecko = GeckoModel(model, biomass_reaction_id='BIOMASS_Ec_iML1515_core_75p37M')
    simulProb = GeckoSimulationProblem(gecko, constraints=constraints)

    simulProb.wt_concentrations = simulProb.simulate().get_protein_concentrations()
    criticalProteins = []
    proteins = [x for x in simulProb.model.proteins if
                x not in criticalProteins and x not in simulProb.objective.keys()]
    proteins.sort()
    LEVELS = [1e-3, 1e-2, 1e-1, 0.5, 0, 1, 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
    decoder = DecoderProtUnderOverExpression(proteins, levels=LEVELS)
    candidate = set(decoder.decode_candidate_ids_to_index(identifiers=UO))
    override = decoder.get_override_simul_problem(candidate=candidate, simulProblem=simulProb)
    res = simulProb.simulate(override)
    print('product', res.ssFluxesDistrib['EX_succ_e'])
    print('biomass', res.ssFluxesDistrib['BIOMASS_Ec_iML1515_core_75p37M'])


def simulate_multi(prot_measure_fractions=None, prot_measure_ggdw=None, constraints={}):
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)



    simulProb = GeckoSimulationProblem(model, constraints=constraints, objective={'r_2056':1}, minimize=True)
    #essencial = simulProb.find_essential_proteins()
    #print(essencial)


    # results = simulProb.simulate()
    # bio = results.ssFluxesDistrib['r_2111']
    # succ_bio = results.ssFluxesDistrib['r_2056']
    # #
    # # print(results.ssFluxesDistrib)
    # #
    # print('biomass', bio)
    # print('product', succ_bio)
    # return None

def define_wt_values(prot_measure_fractions=None, prot_measure_ggdw=None, constraints={}):
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)

    simulProb = GeckoSimulationProblem(model, constraints=constraints)#, method_simul='pFBA')

    wt = simulProb.simulate()
    wt_concentrations = wt.get_protein_concentrations()

    # res = simulProb.simulate()
    # print(res.ssFluxesDistrib['r_2111'], 'biomass')
    # print(res.ssFluxesDistrib['r_2056'], 'succ')


    prot_over = ['P17505', 'P22133', 'P11154', 'P28240', 'P08417', 'P32327', 'P30952', 'P21826', 'P10963', 'P32614', 'P32419']

    for prot in prot_over:
        if wt_concentrations[prot] == 0:
            print(prot)
            wt_concentrations[prot] = 0.00000001

    return wt_concentrations


def get_uniprot_name(prot):
    from Bio import SeqIO
    import urllib.request as url
    handle = url.urlopen('http://www.uniprot.org/uniprot/'+prot+'.xml')
    record = SeqIO.read(handle, "uniprot-xml")
    return record.description

def convert_mmol_to_g(filename):
    from geckopy.data import PROTEIN_PROPERTIES
    df = pandas.read_csv(filename, header=None, index_col=0)
    mmol_series = df.iloc[:, 0]

    grams = []
    for ind, value in mmol_series.iteritems():
        gram = value/1000*PROTEIN_PROPERTIES.loc[ind, 'mw']
        grams.append(gram)

    g_series = pandas.Series(data=grams, index=mmol_series.index)

    return g_series




if __name__ =="__main__":
    import time
    const_aero = {'r_1714_REV':(0,10)}#, 'r_0659_REVNo1':(0,0)}
    const_anaero = {'r_1714_REV':(0,10), 'r_1992_REV':(0,0), 'r_1994_REV':(0,1000), 'r_2189_REV':(0,1000), 'r_1757_REV':(0,1000)}



    const_anaerobic = {'r_1714_REV':(0,10), 'r_1992_REV':(0,0), 'r_1994_REV':(0,1000), 'r_2189_REV':(0,1000), 'r_1757_REV':(0,1000),
                       'r_2137_REV': (0, 1000),'r_2134_REV':(0,1000),'r_1915_REV':(0,1000), 'r_2106_REV': (0,1000),
                       'r_0487No1':(0,0), 'r_0472No1':(0,0),
                       'r_4046':(0,0)}

    ggdw = convert_mmol_to_g('sanchez-mmol_gdw.csv')

    #wt_concentrations = define_wt_values(constraints=const_anaerobic)#, prot_measure_ggdw=ggdw)
    #print(wt.sum())
    #wt_concentrations = wt
    #wt_fluxes = wt[1]

    # print(wt_concentrations['P17505']) #sim
    # print(wt_concentrations['P22133']) #sim
    # print(wt_concentrations['P11154']) #0
    # print(wt_concentrations['P28240']) #0
    # print(wt_concentrations['P08417']) #sim
    # print(wt_concentrations['P32327']) #sim
    # print(wt_concentrations['P30952']) #0
    # print(wt_concentrations['P21826']) #0
    # print(wt_concentrations['P10963']) #0
    # print(wt_concentrations['P32614']) #0
    # print(wt_concentrations['P32419']) #0

    #print(wt_concentrations)

    #simulate_multi(constraints=const_aero, prot_measure_ggdw=ggdw)

    # const_BEST_biomass = {'r_1714_REV':(0,1), 'r_2111':(0.009238,1000),
    #                       'draw_prot_P41939' : (0, 0),
    #                       'draw_prot_P00890' : (2.310152716700675e-06, 10000),
    #                       'draw_prot_P11412' : (0, 0),
    #                       'draw_prot_P37298' : (0, 0),
    #                       'draw_prot_P21954' : (0, 8.621606887004427e-08),
    #                       'draw_prot_P28834' : (0, 0.0)}


    #SUCC
    minBiomassValue_aero = 0.03135
    minBiomassValue_anaero = 0.01875
    minBiomassValue_anaero2 = 0.0169
    minBiomassValue_aero_multipool = 0.00192
    minBiomassValue_anaero_multipool = 0.0009055

    # default_value = (0.00001,10000)
    # const_aero = {'r_1714_REV': (0, 10), 'draw_prot_P32327': default_value,
    #               'draw_prot_P30952': default_value,
    #               'draw_prot_P21826': default_value,
    #               'draw_prot_P10963': default_value,
    #               'draw_prot_P32614': default_value,
    #               'draw_prot_P32419': default_value}

    #KO = [['P33330', 'P00360', 'P06208','P41939', 'Q12680', 'P37303', 'P54115' ]]

    #simulate_KO(target='r_2056', KOs=KO, constraints=const_aero, minBiomassValue=0)
    #print(verify_sol(solution=KO[0], target='r_2056', constraints=const_aero, KO=True, func='wyield', minBiomassValue=minBiomassValue_aero_multipool))

    #uo = [{'P41939':(0,0), 'P00890':(0,10), 'P11412':(0,0), 'P37298':(0,0), 'P21954':(0,0.1), 'P28834':(0,0.5)}]

    # level = 4
    # uo = [{'P12695': (0, 2), 'P28240': (0, 8), 'P41939': (0, 16), 'P00359': (0, 2), 'P11154': (0, 2), 'P20967': (0, 8), 'P16451': (0, 4), 'P46971': (0, 32), 'Q12233': (0, 4), 'P16387': (0, 4), 'Q05584': (0, 16), 'P17423': (0, 0.0625), 'P09624': (0, 4), 'P19262': (0, 2), 'P32614': (0, 4), 'P22580': (0, 0.25), 'P53312': (0, 8), 'P16862': (0, 4)},
    #       {'P46971': (0, 32), 'Q12233': (0, 4), 'P28240': (0, 8), 'P12695': (0, 2), 'P41939': (0, 16), 'P16387': (0, 4), 'Q05584': (0, 16), 'P17423': (0, 0.0625), 'P00359': (0, 2), 'P11154': (0, 2), 'P09624': (0, 4), 'P19262': (0, 2), 'P32614': (0, 4), 'P20967': (0, 8), 'P16451': (0, 4), 'P53312': (0, 8), 'P16862': (0, 4)},
    #       {'P37292': (0, 2), 'P33330': (0, 8), 'Q12233': (0, 4), 'P12695': (0, 2), 'P28240': (0, 8), 'P16387': (0, 4),
    #        'P30952': (0, 16), 'P17423': (0, 0.0625), 'P00359': (0, 2), 'P11154': (0, 2), 'P09624': (0, 4),
    #        'P19262': (0, 2), 'P38113': (0, 2), 'P32614': (0, 4), 'Q12189': (0, 8), 'P20967': (0, 8), 'P16451': (0, 4),
    #        'P16862': (0, 4)},
    #       {'P46971': (0, 32), 'Q12452': (0, 0.125), 'P28240': (0, 8), 'P47164': (0, 0.125), 'Q07804': (0, 0.25),
    #        'P00359': (0, 2), 'P09624': (0, 4), 'P38113': (0, 2), 'P32614': (0, 4), 'P49775': (0, 0.5)}
    #       ]
    #
    #
    # UO = [{
    #     'P17505': (0, level),
    #     'P22133': (0, level),
    #     'P11154': (0, level),
    #     'P28240': (0, level),
    #     'P08417': (0, level),
    #     'P32327': (0, level),
    #     'P30952': (0, level),
    #     'P21826': (0, level),
    #     'P10963': (0, level),
    #     'P32614': (0, level),
    #     'P32419': (0, level),
    #     'P38113': (0,0),
    #     'P00330': (0,0)}]
    #     # 'P26263':(0,0),
    #      #'P16467':(0,0)}]
    #      #'P06169':(0,0)}]


   # UO = [{'Q12122': (0, 2), 'P38858': (0, 4), 'P09624': (0, 4),  'P54115': (0, 0.125), 'P00127': (0, 0.03125), 'P28240': (0, 32), 'P31116': (0, 0.125), 'P43567': (0, 0.5)}]


    # UO = [{'P11412': (0, 0), 'P28834': (0, 0.5), 'P41939': (0, 0), 'P21954': (0, 0.1), 'P37298': (0, 0), 'P00890': (0, 10)}]

    UO = [{'P41939':(0,10)}]

    # print('max biomass')
    bio = simulate_UO(target='r_2056', wt_concentrations= [], wt_fluxes=[], UOs=[{}], constraints=const_aero, minBiomassValue=minBiomassValue_aero)#, prot_measure_ggdw=ggdw)

    # print('99%')
    # const_anaerobic['r_2111'] = (bio*0.99, 1000)
    #
    # print('max succinate')
    # simulate_UO(target='r_2056', wt_concentrations= [], wt_fluxes=[], UOs=UO, constraints=const_anaerobic,
    #              objective={'r_2056':1})#, prot_measure_ggdw=ggdw)
    # print('min succinate')
    # simulate_UO(target='r_2056',  wt_concentrations= wt_concentrations, wt_fluxes=[], UOs=UO, constraints=const_anaerobic,
    #              minBiomassValue=minBiomassValue_anaero2, objective={'r_2056':1}, minimize=True)#, prot_measure_ggdw=ggdw)
    #
    # print('90%')
    # const_anaerobic['r_2111'] = (bio *0.90, 10000)
    # print('max succinate')
    # simulate_UO(target='r_2056', wt_concentrations=wt_concentrations, wt_fluxes=wt_fluxes, UOs=UO, constraints=const_anaerobic,
    #              minBiomassValue=minBiomassValue_aero, objective={'r_2056': 1}, prot_measure_ggdw=ggdw)
    # print('min succinate')
    # simulate_UO(target='r_2056', wt_concentrations=wt_concentrations, wt_fluxes=wt_fluxes, UOs=UO, constraints=const_anaerobic,
    #              minBiomassValue=minBiomassValue_aero, objective={'r_2056': 1}, minimize=True, prot_measure_ggdw=ggdw)
    #
    # print('10%')
    # const_anaerobic['r_2111'] = (bio *0.10, 10000)
    # print('max succinate')
    # simulate_UO(target='r_2056', wt_concentrations=wt_concentrations, wt_fluxes=[], UOs=UO, constraints=const_anaerobic,
    #              minBiomassValue=minBiomassValue_aero, objective={'r_2056': 1})# prot_measure_ggdw=ggdw)
    #
    # print('min succinate')
    # simulate_UO(target='r_2056', wt_concentrations=wt_concentrations, wt_fluxes=[], UOs=UO, constraints=const_anaerobic,
    #              minBiomassValue=minBiomassValue_aero, objective={'r_2056': 1}, minimize=True)#, prot_measure_ggdw=ggdw)

    # print('UO')
    # simulate_UO(target='r_2056', wt_concentrations=wt_concentrations, UOs=uo, constraints=const_aero, minBiomassValue=minBiomassValue_aero_multipool, prot_measure_ggdw=ggdw)
    # print(verify_sol(solution=uo, target='r_2056', constraints=const_aero, KO=False, func='wyield', prot_measure_ggdw=ggdw, minBiomassValue=minBiomassValue_aero_multipool))

    # # TYR

    # start = time.time()
    # simulate_KO(target='r_1913', KO=listas, constraints=const)
    # #verify_sol(candidate=lista, target='r_1913', constraints=const, KO=True, func='TF')
    # end = time.time()
    # print(end-start, ' seconds')

    # # PHE

    # start = time.time()
    # simulate_KO(target='r_1903', KO=listas, constraints=const)
    # end = time.time()
    # print(end-start, ' seconds')
    # verify_sol(candidate=lista, target='r_1903', constraints=const, KO=True, func='TF')

    # # NEW MODEL - ECOLI
    # lista = ['P52043', 'P0AC44', 'P0AGE9', 'P0A7E1']
    # dic = {'P00864':(0,10)}
    # const_ecoli = {'EX_glc__D_e_REV': (0, 1)}
    # ecoli_model_KO(KO=[], constraints=const_ecoli)
    # ecoli_model_UO(UO=dic, constraints=const_ecoli)
