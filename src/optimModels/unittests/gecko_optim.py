from geckopy import GeckoModel
from optimModels.optimization.evaluation_functions import build_evaluation_function
from optimModels.simulation.simul_problems import GeckoSimulationProblem
from optimModels.optimization.run import gecko_strain_optim
from optimModels.utils.constantes import optimType
from cobra.io import read_sbml_model

basePath = "C:/Users/BiSBII/Results/"
#LEVELS = [1e-3, 1e-2, 1e-1, 0.5, 0, 1 , 5, 10, 50, 1e2, 5e2, 1e3, 1e4]
LEVELS = [0, 2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5]

# Optimization:
# model: GECKO YEAST
# Obj Func: TargetFlux
# type: Protein Knockouts
# product: succinate

def prot_ko_optim_succ_TF (fileRes, prot_measure_fractions =None, prot_measure_ggdw= None, constraints = None, isMultiProc=False, size=5, iniPop = None):

    #load model
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)


    simulProb = GeckoSimulationProblem(model, constraints= constraints, method_simul='pFBA')

    #essenciais = simulProb.find_essential_proteins() #this function takes a lot of time to run (6000 proteins)
    #print(essenciais)
    essenciais = ['P04076', 'P15624', 'P14020', 'P40545', 'P38708', 'P50861', 'P39692', 'P39954', 'P46655', 'P04173',
                  'P47169', 'P05375', 'P03965', 'P06197',
                  'P31116', 'P16120', 'P38707', 'P00815', 'P07806', 'P07702', 'P18544', 'Q05911', 'P15625', 'Q12055',
                  'P32347', 'Q01217', 'P33312', 'P07149',
                  'P38221', 'P07264', 'Q02196', 'P29704', 'P22768', 'P40825', 'P41940', 'P20051', 'P07807', 'P54839',
                  'P07285', 'P15180', 'Q05506', 'P08465',
                  'P38066', 'P36148', 'P53045', 'Q99258', 'Q12362', 'P08524', 'P15700', 'P32288', 'P28777', 'P04802',
                  'P19414', 'P28272', 'P53199', 'P13188',
                  'P11353', 'P08456', 'P06168', 'P41338', 'P25340', 'P07283', 'P38999', 'P53852', 'P07277', 'P05373',
                  'P16622', 'P09950', 'P33734', 'P00899',
                  'P07263', 'Q04728', 'P38145', 'P10614', 'P00958', 'P05694', 'P08566', 'P32263', 'P38891', 'P29952',
                  'P15496', 'P32476', 'P39522', 'P16603',
                  'P07172', 'P48445', 'P17423', 'P38635', 'P50113', 'P49367', 'P24521', 'P38998', 'P32462', 'P07342',
                  'P00912', 'P07258', 'P47176', 'P05150',
                  'Q12122', 'P06633', 'P07244', 'P40495', 'P21147', 'P18408', 'P07259', 'P36421', 'P38604', 'P19097',
                  'P10869', 'P04801', 'P32178', 'P09436',
                  'Q12452', 'P00498', 'P15454', 'P06106', 'P00937', 'P00931', 'Q12109', 'P38972', 'P06785', 'P38625',
                  'P13663', 'P26637', 'Q12189', 'P03962',
                  'P31688', 'P11986', 'Q00955', 'P20049', 'Q00764', 'P32377', 'P28789', 'P29509', 'P32452']

    evalFunc = build_evaluation_function("targetFlux", ["r_2056"])

    gecko_strain_optim(simulProblem=simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc, candidateSize= size,
                       resultFile=fileRes, criticalProteins=essenciais, initPopFile=iniPop) #KO_Reaction by default


# Optimization:
# model: GECKO YEAST
# Obj Func: BPCY
# type: Protein Knockouts
# product: succinate

def prot_ko_optim_succ_BPCY(fileRes, prot_measure_fractions=None, prot_measure_ggdw=None, constraints=None, isMultiProc=False, size=5, iniPop=None):
    # load model
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)

    simulProb = GeckoSimulationProblem(model, constraints=constraints)

    essenciais = simulProb.find_essential_proteins()  # this function takes a lot of time to run (6000 proteins)

    evalFunc = build_evaluation_function("BPCY", "r_2111", "r_2056", "r_1714_REV")  # max succ exchange

    gecko_strain_optim(simulProblem=simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc,
                       candidateSize=size, resultFile=fileRes, criticalProteins=essenciais, initPopFile=iniPop)  # KO_Reaction by default


# Optimization:
# model: GECKO YEAST
# Obj Func: WYIELD
# type: Protein Knockouts
# product: succinate

def prot_ko_optim_succ_WYIELD (fileRes, prot_measure_fractions =None, prot_measure_ggdw= None, constraints = None, isMultiProc=False, size=5, iniPop = None):

    #load model
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)


    simulProb = GeckoSimulationProblem(model, constraints= constraints)


    #essenciais = simulProb.find_essential_proteins() #this function takes a lot of time to run (6000 proteins)
    #print(essenciais)
    '''
    essenciais = ['P04076', 'P15624', 'P14020', 'P40545', 'P38708', 'P50861', 'P39692', 'P39954', 'P46655', 'P04173',
                  'P47169', 'P05375', 'P03965', 'P06197',
                  'P31116', 'P16120', 'P38707', 'P00815', 'P07806', 'P07702', 'P18544', 'Q05911', 'P15625', 'Q12055',
                  'P32347', 'Q01217', 'P33312', 'P07149',
                  'P38221', 'P07264', 'Q02196', 'P29704', 'P22768', 'P40825', 'P41940', 'P20051', 'P07807', 'P54839',
                  'P07285', 'P15180', 'Q05506', 'P08465',
                  'P38066', 'P36148', 'P53045', 'Q99258', 'Q12362', 'P08524', 'P15700', 'P32288', 'P28777', 'P04802',
                  'P19414', 'P28272', 'P53199', 'P13188',
                  'P11353', 'P08456', 'P06168', 'P41338', 'P25340', 'P07283', 'P38999', 'P53852', 'P07277', 'P05373',
                  'P16622', 'P09950', 'P33734', 'P00899',
                  'P07263', 'Q04728', 'P38145', 'P10614', 'P00958', 'P05694', 'P08566', 'P32263', 'P38891', 'P29952',
                  'P15496', 'P32476', 'P39522', 'P16603',
                  'P07172', 'P48445', 'P17423', 'P38635', 'P50113', 'P49367', 'P24521', 'P38998', 'P32462', 'P07342',
                  'P00912', 'P07258', 'P47176', 'P05150',
                  'Q12122', 'P06633', 'P07244', 'P40495', 'P21147', 'P18408', 'P07259', 'P36421', 'P38604', 'P19097',
                  'P10869', 'P04801', 'P32178', 'P09436',
                  'Q12452', 'P00498', 'P15454', 'P06106', 'P00937', 'P00931', 'Q12109', 'P38972', 'P06785', 'P38625',
                  'P13663', 'P26637', 'Q12189', 'P03962',
                  'P31688', 'P11986', 'Q00955', 'P20049', 'Q00764', 'P32377', 'P28789', 'P29509', 'P32452']
    '''
    essenciais = ['P11986', 'P28777', 'P06785', 'P15700', 'P15624', 'P00899', 'P22768', 'Q01217', 'P29509', 'P10614', 'P47176', 'P46655', 'P32263', 'P07806', 'P38891', 'P39692', 'P28272', 'Q05506', 'P07342', 'P16120', 'Q00764', 'P07149', 'Q12122', 'P07283', 'Q00955', 'P33734', 'Q99258', 'P00958', 'P40825', 'P15180', 'P20051', 'P36148', 'P32178', 'Q12362', 'P08456', 'P32476', 'P07264', 'P13663', 'P07258', 'P38145', 'P29704', 'P38708', 'P10869', 'P17423', 'P05375', 'P40495', 'P18544', 'P05150', 'P00912', 'P15496', 'P41338', 'P29952', 'P47169', 'P39522', 'P06168', 'P50861', 'P38604', 'P13188', 'P49367', 'P00498', 'P15625', 'Q12452', 'P32377', 'P09436', 'P48445', 'P03965', 'P08465', 'P32452', 'P14020', 'P00937', 'P21147', 'P53045', 'P04076', 'P53852', 'P07807', 'P25340', 'P31688', 'P19097', 'P36421', 'P06106', 'P07172', 'P06633', 'P00931', 'P26637', 'P07259', 'P08566', 'P54839', 'P38998', 'P38635', 'P20049', 'P07285', 'P19414', 'P07263', 'P16603', 'P04802', 'P39954', 'P07277', 'P06197', 'P38999', 'Q12109', 'P24521', 'Q02196', 'P00815', 'P33312', 'Q04728', 'P38707', 'P07702', 'P04801', 'P32462', 'P05694', 'P31116', 'P41940', 'P38221', 'P08524', 'P40545', 'P15454', 'P50113', 'P04173', 'P03962', 'P04161']

    evalFunc = build_evaluation_function("WYIELD", "r_2111", "r_2056", alpha=0.3, minBiomassValue=0.00192)

    gecko_strain_optim(simulProblem=simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc, candidateSize= size,
                        resultFile=fileRes, criticalProteins=essenciais, initPopFile=iniPop) #KO_Reaction by default


# Optimization:
# model: GECKO YEAST
# Obj Func: TargetFlux
# type: Protein UnderOverExpression
# product: succinate

def prot_uo_optim_succ_TF (fileRes, prot_measure_fractions =None, prot_measure_ggdw= None,  constraints = None, isMultiProc=False, size=5, iniPop = None):
    #load model
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)

    simulProb = GeckoSimulationProblem(model=model, constraints=constraints, minimize=True, objective={'r_2056':1})
    wt_fluxes = simulProb.simulate()
    simulProb.wt_concentrations = wt_fluxes.get_protein_concentrations()

    proteins = ['P38113', 'P00330', 'P00331', 'P10127', 'P07246', 'Q12358', 'P38715', 'P32327', 'P10963', 'P17505',
                'P08417', 'P32614', 'P28240', 'P21826', 'P11154', 'P22133', 'P30952', 'P32419', 'P32191', 'Q00055',
                'P41911', 'P26263', 'P07342', 'P37298', 'P47052', 'P33421', 'P21801', 'Q00711', 'P28834', 'P28241',
                'P41939', 'P53982', 'P21954', 'P09624', 'P19262', 'P20967', 'P53312', 'P53598', 'P26263', 'P06169',
                'P16467', 'P00890', 'P43635', 'P08679', 'P11412', 'P42941', 'P33330', 'P40054', 'P40510', 'P32771',
                'P25377', 'Q12458', 'Q04894']

    all_prot = simulProb.model.proteins
    essencial = set(all_prot) - set(proteins)


    evalFunc = build_evaluation_function("targetFlux", ["r_2056"])

    gecko_strain_optim(simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc, candidateSize=size,
                       resultFile=fileRes, type=optimType.PROTEIN_UO, criticalProteins=essencial, initPopFile=iniPop)


# Optimization:
# model: GECKO YEAST
# Obj Func: BPCY
# type: Protein UnderOverExpression
# product: succinate

def prot_uo_optim_succ_BPCY(fileRes, prot_measure_fractions=None, prot_measure_ggdw=None, constraints=None,
                          isMultiProc=False, size=5, iniPop=None):
    # load model
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)


    simulProb = GeckoSimulationProblem(model, constraints=constraints)

    evalFunc = build_evaluation_function("BPCY", "r_2111", "r_2056", "r_1714_REV")  # max succ exchange

    gecko_strain_optim(simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc, candidateSize=size,
                       resultFile=fileRes, type=optimType.PROTEIN_UO, initPopFile=iniPop)

# Optimization:
# model: GECKO YEAST
# Obj Func: WYIELD
# type: Protein UnderOverExpression
# product: succinate

def prot_uo_optim_succ_WYIELD(fileRes, prot_measure_fractions=None, prot_measure_ggdw=None, constraints=None,
                          isMultiProc=False, size=5, iniPop=None):
    # load model
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)


    simulProb = GeckoSimulationProblem(model, constraints=constraints)
    wt_fluxes = simulProb.simulate()
    simulProb.wt_concentrations = wt_fluxes.get_protein_concentrations()
    simulProb._protein_rev_reactions = simulProb.get_protein_rev_reactions()
    #print(simulProb.wt_concentrations)



    evalFunc = build_evaluation_function("WYIELD", "r_2111", "r_2056", alpha=0.3, minBiomassValue=0.0001)

    gecko_strain_optim(simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc, candidateSize=size,
                        resultFile=fileRes, type=optimType.PROTEIN_UO, initPopFile=iniPop)


# Optimization:
# model: GECKO YEAST
# Obj Func: Target Flux
# type: Protein Knockouts
# product: tyrosine

def prot_ko_optim_tyr_TF (fileRes, prot_measure_fractions =None, prot_measure_ggdw= None, constraints = None, isMultiProc=False, size=5, iniPop = None):

    #load model
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)


    simulProb = GeckoSimulationProblem(model, constraints= constraints)

    #essenciais = simulProb.find_essential_proteins() #this function takes a lot of time to run (6000 proteins)
    essenciais =  ['P04076', 'P15624', 'P14020', 'P40545', 'P38708', 'P50861', 'P39692', 'P39954', 'P46655', 'P04173',
                  'P47169', 'P05375', 'P03965', 'P06197',
                  'P31116', 'P16120', 'P38707', 'P00815', 'P07806', 'P07702', 'P18544', 'Q05911', 'P15625', 'Q12055',
                  'P32347', 'Q01217', 'P33312', 'P07149',
                  'P38221', 'P07264', 'Q02196', 'P29704', 'P22768', 'P40825', 'P41940', 'P20051', 'P07807', 'P54839',
                  'P07285', 'P15180', 'Q05506', 'P08465',
                  'P38066', 'P36148', 'P53045', 'Q99258', 'Q12362', 'P08524', 'P15700', 'P32288', 'P28777', 'P04802',
                  'P19414', 'P28272', 'P53199', 'P13188',
                  'P11353', 'P08456', 'P06168', 'P41338', 'P25340', 'P07283', 'P38999', 'P53852', 'P07277', 'P05373',
                  'P16622', 'P09950', 'P33734', 'P00899',
                  'P07263', 'Q04728', 'P38145', 'P10614', 'P00958', 'P05694', 'P08566', 'P32263', 'P38891', 'P29952',
                  'P15496', 'P32476', 'P39522', 'P16603',
                  'P07172', 'P48445', 'P17423', 'P38635', 'P50113', 'P49367', 'P24521', 'P38998', 'P32462', 'P07342',
                  'P00912', 'P07258', 'P47176', 'P05150',
                  'Q12122', 'P06633', 'P07244', 'P40495', 'P21147', 'P18408', 'P07259', 'P36421', 'P38604', 'P19097',
                  'P10869', 'P04801', 'P32178', 'P09436',
                  'Q12452', 'P00498', 'P15454', 'P06106', 'P00937', 'P00931', 'Q12109', 'P38972', 'P06785', 'P38625',
                  'P13663', 'P26637', 'Q12189', 'P03962',
                  'P31688', 'P11986', 'Q00955', 'P20049', 'Q00764', 'P32377', 'P28789', 'P29509', 'P32452']

    evalFunc = build_evaluation_function("targetFlux", ["r_1913"])

    gecko_strain_optim(simulProblem=simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc, candidateSize= size,
                       resultFile=fileRes, criticalProteins=essenciais, initPopFile=iniPop) #KO_Reaction by default

# Optimization:
# model: GECKO YEAST
# Obj Func: Target Flux
# type: Protein Knockouts
# product: phenylalanine

def prot_ko_optim_phe_TF (fileRes, prot_measure_fractions =None, prot_measure_ggdw= None, constraints = None, isMultiProc=False, size=5, iniPop = None):

    #load model
    if prot_measure_fractions is None and prot_measure_ggdw is None:
        model = GeckoModel("single-pool")
    else:
        model = GeckoModel("multi-pool")
        if prot_measure_fractions:
            model.limit_proteins(fractions=prot_measure_fractions)
        else:
            model.limit_proteins(ggdw=prot_measure_ggdw)


    simulProb = GeckoSimulationProblem(model, constraints= constraints)

    essenciais = simulProb.find_essential_proteins() #this function takes a lot of time to run (6000 proteins)

    evalFunc = build_evaluation_function("targetFlux", ["r_1903"])

    gecko_strain_optim(simulProblem=simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc, candidateSize= size,
                        resultFile=fileRes, criticalProteins=essenciais, initPopFile=iniPop) #KO_Reaction by default

# Optimization:
# model: GECKO ECOLI
# Obj Func: Target Flux
# type: Protein Knockouts
# product: succinate

def ecoli_ko_succ (fileRes, prot_measure_fractions =None, prot_measure_ggdw= None, constraints = None, isMultiProc=False, size=5, iniPop = None):
    SBML_FILE = '../../../examples/models/eciML1515_batch.xml'
    cobra = read_sbml_model(SBML_FILE)

    model = GeckoModel(cobra, biomass_reaction_id='BIOMASS_Ec_iML1515_core_75p37M')

    simulProb = GeckoSimulationProblem(model, constraints=constraints)


    #essenciais = simulProb.find_essential_proteins()  # this function takes a lot of time to run (6000 proteins)
    #print(essenciais)
    essenciais = ['P00904', 'P0A6T5', 'P60752', 'P0A790', 'P00547', 'P06986', 'P12008', 'Q46893', 'P0AGK1', 'P0A786', 'P11447', 'P62617', 'P00934', 'P00909', 'P00903', 'P0A8K1', 'P04951', 'P0AB89', 'P17443', 'P0ABZ4', 'P15254', 'P0A9J8', 'P26647', 'P0A6I3', 'P0AGG0', 'P0A725', 'P45578', 'P0A6D3', 'P0A6X1', 'P0AC75', 'P0A877', 'P0A6W3', 'P17952', 'P0A6K1', 'P0ABD8', 'P0A7B5', 'P07004', 'P0A6C8', 'P0ABD5', 'P09126', 'P23908', 'P0A9V1', 'P05020', 'P62623', 'P05194', 'P19624', 'P04036', 'P22634', 'P05041', 'P0AD68', 'P0AFU8', 'P0ABQ0', 'P0AB80', 'P0A720', 'P0A6L2', 'P0A749', 'P29680', 'P0AF16', 'P0ABH7', 'P60757', 'P0A752', 'P23830', 'P60472', 'P43341', 'P17854', 'P06989', 'P0A722', 'P0A6C5', 'P0AF98', 'P30126', 'P0ADV1', 'P60595', 'P22188', 'P31554', 'P06987', 'P0A6A6', 'P60546', 'P0A6I6', 'P0A9Q5', 'P0AEA8', 'P05793', 'P0A9D4', 'P22939', 'P15770', 'P0A9D8', 'P07623', 'P24182', 'P00935', 'P0A7I7', 'P15639', 'P0AC16', 'P26281', 'P05791', 'P0AD57', 'P07023', 'P0ADV9', 'P0ACB4', 'P10902', 'P14900', 'P62615', 'P0ABG4', 'P10441', 'P60664', 'P08373', 'P06983', 'P0AC13', 'P0A7E3', 'P0A6E4', 'P61714', 'P0AF12', 'P10371', 'P31120', 'P00895', 'P0A884', 'P0ADC6', 'P0ADC1', 'P0A6A8', 'P08192', 'P04805', 'P0ABG1', 'P21645', 'P25539', 'P0AG40', 'P0A6J1', 'P00861', 'P0A715', 'P08200', 'P0A7J0', 'P30011', 'P06988', 'P0A7E5', 'P08244', 'P0A879', 'P0AAI9', 'P23845', 'P0AD65', 'P77488', 'P0A794', 'P62620', 'P17169', 'P0AEK2', 'P11446', 'P27300', 'P0A9L8', 'P0A9Q9', 'P04079', 'P11880', 'P0AED7', 'P0ACC7', 'P23893', 'P0AEZ1', 'P0ACB2', 'P0A817', 'P09151', 'P31663', 'P45568']

    evalFunc = build_evaluation_function("targetFlux", ["EX_succ_e"])

    gecko_strain_optim(simulProblem=simulProb, evaluationFunc=evalFunc, levels=LEVELS, isMultiProc=isMultiProc,
                       candidateSize=size, resultFile=fileRes, criticalProteins=essenciais, initPopFile=iniPop)  # KO_Reaction by default


def convert_mmol_to_g(filename):
    import pandas
    from geckopy.data import PROTEIN_PROPERTIES
    df = pandas.read_csv(filename, header=None, index_col=0)
    mmol_series = df.iloc[:, 0]

    grams = []
    for ind, value in mmol_series.iteritems():
        gram = value/1000*PROTEIN_PROPERTIES.loc[ind, 'mw']
        grams.append(gram)

    g_series = pandas.Series(data=grams, index=mmol_series.index)

    return g_series




if __name__ == '__main__':
    import time
    ############################################# YEAST ##############################################################

    const_TF = {'r_1714_REV': (0, 10)}
    const_anaerobic = {'r_1714_REV': (0, 10), 'r_1992_REV': (0, 0), 'r_1994_REV': (0, 1000), 'r_2189_REV': (0, 1000),
                       'r_1757_REV': (0, 1000)}

    # const_BPCY = {'r_1714_REV': (0, 1)}
    size = 10
    ggdw = convert_mmol_to_g('sanchez-mmol_gdw.csv')

    # #SUCC_KO_TF
    fileRes = basePath + "multi-pool_test_ko_ana.csv"
    start = time.time()
    prot_ko_optim_succ_WYIELD(fileRes=fileRes, prot_measure_ggdw=ggdw,constraints=const_anaerobic, size=size)
    #prot_uo_optim_succ_WYIELD(fileRes=fileRes, prot_measure_ggdw=ggdw,constraints=const_anaerobic, size=size)
    #prot_uo_optim_succ_WYIELD(fileRes=fileRes, constraints=const_TF, size=size)
    end = time.time()
    #print('it took ', end - start, ' seconds')

    # #SUCC_KO_BPCY
    # fileRes = basePath + "Gecko_KO_SUCC_size5_BPCY.csv"
    # prot_ko_optim_succ_BPCY(fileRes=fileRes, constraints=const_BPCY, size=size)

    #SUCC_UO_TF
    # start = time.time()
    # fileRes = basePath + "WYIELD_Gecko_UO_SUCC_size8.csv"
    # lista = [10]
    # for x in lista:
    #     size = x
    #     fileRes = basePath + 'scale_WYIELD_Gecko_UO_SUCC_size'+str(size)+'.csv'
    #     prot_uo_optim_succ_WYIELD(fileRes=fileRes, constraints=const_TF, size=size)
    # # # end = time.time()
    # print('it took ',end-start,' seconds')

    # #SUCC_UO_BPCY
    # fileRes = basePath + "Gecko_UO_SUCC_size5_BPCY.csv"
    # prot_uo_optim_succ_BPCY(fileRes=fileRes, constraints=const_BPCY, size=size)
    #
    # #TYR_KO_TF
    #fileRes = basePath + "pFBA_Gecko_KO_TYR_size5_TF.csv"
    #size=6
    #prot_ko_optim_tyr_TF(fileRes=fileRes, constraints=const_TF, size=size)


    ############################################# ECOLI ##############################################################
    # const_TF = {'EX_glc__D_e_REV': (0, 1), 'EX_o2_e_REV'}
    # lista = [4,5,6,7,8]
    # for x in lista:
    #      size = x
    #      fileRes = basePath + "GeckoEcoli_KO_SUCC_"+str(size)+"_TF.csv"
    #      ecoli_ko_succ(fileRes=fileRes, constraints=const_TF, size=size)