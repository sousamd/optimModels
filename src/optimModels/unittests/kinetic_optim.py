from optimModels.optimization.evaluation_functions import build_evaluation_function
from optimModels.optimization.run import kinetic_strain_optim
from optimModels.model.kineticModel import load_kinetic_model
from optimModels.simulation.simul_problems import KineticSimulationProblem
from optimModels.utils.configurations import KineticConfigurations
from optimModels.optimization.evolutionary_computation import EAConfigurations
from collections import OrderedDict
from optimModels.utils.constantes import optimType
import warnings


LEVELS = [0, 2 ** -5, 2 ** -4, 2 ** -3, 2 ** -2, 2 ** -1, 2 ** 1, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5]
#basePath = "/home/scorreia/Decaf/"
basePath = "C:/Users/BISBII/Results/"

def ko_chassagnole(isMultiProc=False, size=1):
    sbmlFile = '../../../examples/models/chassagnole2002.xml'
    fileRes = basePath + "optim_Chassagnole_Serine_ko2_"+str(size)+".csv"


    model1 = load_kinetic_model(sbmlFile, [])
    #print(model.parsedRates)

    parameters = list(model1.get_parameters())
    params = parameters[10:]

    map = {}
    for param in params:
        reaction = param.split('_')[0]
        if reaction not in map.keys():
            map[reaction] = [param]

    print(map)

    model = load_kinetic_model(sbmlFile, map)


    objFunc = build_evaluation_function("targetFlux", ["vsersynth"])
    simulProblem = KineticSimulationProblem(model, tSteps=[0, KineticConfigurations.STEADY_STATE_TIME])
    result = kinetic_strain_optim(simulProblem, objFunc=objFunc, levels=None, type= optimType.REACTION_KO, isMultiProc=True, candidateSize=size, resultFile=fileRes)
    result.print()


    #SIMULATION
    '''
    over = {'vSynth1_rmaxSynth1':0}
    ko = OrderedDict(over)
    from optimModels.simulation.override_simul_problem import OverrideKineticSimulProblem
    override = OverrideKineticSimulProblem(ko)
    res = simulProblem.simulate(overrideSimulProblem=override)
    res.print()
    '''


def under_over_chassagnole(isMultiProc=False, size=1):
    sbmlFile = '../../../examples/models/chassagnole2002.xml'
    fileRes = basePath + "optim_Chassagnole_Serine_underover2.csv"

    model = load_kinetic_model(sbmlFile,[])

    parameters = list(model.get_parameters())
    params = parameters[10:]

    map = {}
    for param in params:
        reaction = param.split('_')[0]
        if reaction not in map.keys():
            map[reaction] = [param]

    model = load_kinetic_model(sbmlFile, map)

    objFunc = build_evaluation_function("targetFlux", ["vsersynth"])
    simulProblem = KineticSimulationProblem(model, tSteps=[0, KineticConfigurations.STEADY_STATE_TIME])
    result = kinetic_strain_optim(simulProblem, objFunc=objFunc, levels=LEVELS, criticalParameters=[], isMultiProc=False, candidateSize=size, resultFile=fileRes, type=optimType.REACTION_UO)
    result.print()


def ko_jahan(isMultiProc=False, size=1):
    sbmlFile = '../../../examples/models/Jahan2016_chemostat_fixed.xml'
    fileRes = basePath + "optim_Jahan_Suc_ko.csv"
    mapParamReacs = {"vE_6PGDH": ["v6PGDH_max"], "vE_Ack": ["vAck_max"], "vE_Ack_medium": ["vAck_max"],
                     "vE_Cya": ["vCya_max"], "vE_Eda": ["vEda_max"], "vE_Edd": ["vEdd_max"], "vE_Fum": ["Fum"],
                     "vE_G6PDH": ["vG6PDH_max"], "vE_MDH": ["MDH"], "vE_Pgi": ["vPgi_max"],
                     "vE_Pgl": ["vPgl_max"], "vE_Pta": ["vPta_max"], "vE_R5PI": ["vR5PI_max"], "vE_Ru5P": ["vRu5P_max"],
                     "vE_Tal": ["vTal_max"], "vE_TktA": ["vTktA_max"], "vE_TktB": ["vTktB_max"],
                     "vE_cAMPdegr": ["vcAMPdegr_max"], "vNonPTS": ["vNonPTS_max"], "vNonPTS_medium": ["vNonPTS_max"],
                     "vPTS4": ["vPTS4_max"], "vPTS4_medium": ["vPTS4_max"], "vE_AceKki": ["AceK"],
                     "vE_AceKph": ["AceK"], "vE_Acs": ["Acs"], "vE_Acs_medium": ["Acs"], "vE_CS": ["CS"],
                     "vE_Fba": ["Fba"], "vE_Fbp": ["Fbp"], "vE_GAPDH": ["GAPDH"], "vE_Glk": ["Glk"],
                     "vE_ICDH": ["ICDH"], "vE_Icl": ["Icl"], "vE_MS": ["MS"], "vE_Mez": ["Mez"], "vE_PDH": ["PDH"],
                     "vE_Pck": ["Pck"], "vE_Pfk": ["Pfk"], "vE_Ppc": ["Ppc"], "vE_Pps": ["Pps"], "vE_Pyk": ["Pyk"],
                     "vE_SDH": ["SDH"], "vE_aKGDH": ["aKGDH"]}

    model = load_kinetic_model(sbmlFile, mapParamReacs)

    objFunc = build_evaluation_function("targetFlux", ["vD_SUC"])
    simulProblem= KineticSimulationProblem(model)
    result = kinetic_strain_optim(simulProblem, objFunc=objFunc, levels=LEVELS, isMultiProc=isMultiProc, candidateSize=size, resultFile=fileRes)
    result.print()



def ko_millard(isMultiProc=False, size=1):
    EAConfigurations.MAX_CANDIDATE_SIZE = size;


    sbmlFile = '../../../examples/models/Millard2016v2.xml'
    fileRes = basePath  + "Results/optim_Millard_acet_ko_"+str(size)+".csv"
    #fileLastGen = basePath + "Results/optim_Millard_acet_ko_" + str(size) + "_lastgen.csv"

    mapParamReacs = OrderedDict([('PTS_4', ['eiicbP']), ('PTS_0', ['ei']), ('PTS_1', ['eiP']), ('PTS_2', ['eiia']), ('PTS_3', ['eiicb']),
                      ('PGI', ['PGI_Vmax']), ('PFK', ['PFK_Vmax']), ('FBA', ['FBA_Vmax']), ('TPI', ['TPI_Vmax']),
                      ('GDH', ['GDH_Vmax']), ('PGK', ['PGK_Vmax']), ('GPM', ['GPM_Vmax']), ('ENO', ['ENO_Vmax']),
                      ('PYK', ['PYK_Vmax']), ('ZWF', ['ZWF_Vmax']), ('PGL', ['PGL_Vmax']), ('GND', ['GND_Vmax']),
                      ('RPE', ['RPE_Vmax']), ('RPI', ['RPI_Vmax']), ('X5P_GAP_TKT', ['tkt']), ('F6P_E4P_TKT', ['tktC2']),
                      ('S7P_R5P_TKT', ['tktC2']), ('F6P_GAP_TAL', ['talC3']), ('S7P_E4P_TAL', ['tal']), ('FBP', ['FBP_Vmax']),
                      ('PPC', ['PPC_Vmax']), ('PCK', ['PCK_Vmax']), ('PPS', ['PPS_Vmax']), ('MAD', ['MAD_Vmax']),
                      ('PDH', ['PDH_Vmax']), ('GLT', ['GLT_Vmax']), ('ACN_1', ['ACN_1_Vmax']), ('ACN_2', ['ACN_2_Vmax']),
                      ('ICD', ['icd']), ('LPD', ['LPD_Vmax']), ('SK', ['SK_Vmax']), ('SDH', ['SDH_Vmax']), ('FUMA', ['FUMA_Vmax']),
                      ('MQO', ['MQO_Vmax']), ('MDH', ['MDH_Vmax']), ('ACEA', ['ACEA_Vmax']), ('ACEB', ['ACEB_Vmax']),
                      ('EDD', ['EDD_Vmax']), ('EDA', ['EDA_Vmax']), ('NADH_req', ['NADH_req_Vmax']), ('ATP_syn', ['ATP_syn_Vmax']),
                      ('ACK', ['ACK_Vmax']), ('ACS', ['ACS_Vmax']), ('PTA', ['PTA_Vmax']), ('MYTBO', ['MYTBO_Vmax']),
                      ('SQR', ['SQR_Vmax']), ('NDHII', ['NDHII_Vmax']), ('GROWTH', ['GROWTH_Vmax']), ('ATP_MAINTENANCE', ['ATP_MAINTENANCE_Vmax']),
                      ('XCH_GLC', ['XCH_GLC_Vmax']), ('PIT', ['PIT_Vmax']), ('XCH_P', ['XCH_P_Vmax']), ('XCH_ACE1', ['XCH_ACE1_Vmax']),
                      ('XCH_ACE2', ['XCH_ACE2_Vmax'])])

    model = load_kinetic_model(sbmlFile, mapParamReacs)

    objFunc = build_evaluation_function("targetFlux", ["_ACE_OUT"])
    simulProblem= KineticSimulationProblem(model, tSteps=[0, KineticConfigurations.STEADY_STATE_TIME])

    result = kinetic_strain_optim(simulProblem, objFunc=objFunc, levels=None,
                                  criticalParameters=['ATP_MAINTENANCE_Vmax', 'GROWTH_Vmax', 'NDHII_Vmax', 'PIT_Vmax', 'eiicbP', 'ei', 'eiP', 'eiia'],
                                  isMultiProc=isMultiProc, resultFile=fileRes, initPopFile=None)
    result.print()


def under_over_millard(isMultiProc=False, size=1):
    EAConfigurations.MAX_CANDIDATE_SIZE = size;

    sbmlFile = '../../../examples/models/Millard2016v2.xml'
    fileRes = basePath + "Results/optim_Millard_acet_underover_"+str(size)+".csv"
    fileLastGen = basePath + "Results/optim_Millard_acet_underover_" + str(size) + "_lastgen.csv"
    mapParamReacs = OrderedDict([('PTS_4', ['eiicbP']), ('PTS_0', ['ei']), ('PTS_1', ['eiP']), ('PTS_2', ['eiia']), ('PTS_3', ['eiicb']),
                      ('PGI', ['PGI_Vmax']), ('PFK', ['PFK_Vmax']), ('FBA', ['FBA_Vmax']), ('TPI', ['TPI_Vmax']),
                      ('GDH', ['GDH_Vmax']), ('PGK', ['PGK_Vmax']), ('GPM', ['GPM_Vmax']), ('ENO', ['ENO_Vmax']),
                      ('PYK', ['PYK_Vmax']), ('ZWF', ['ZWF_Vmax']), ('PGL', ['PGL_Vmax']), ('GND', ['GND_Vmax']),
                      ('RPE', ['RPE_Vmax']), ('RPI', ['RPI_Vmax']), ('X5P_GAP_TKT', ['tkt']), ('F6P_E4P_TKT', ['tktC2']),
                      ('S7P_R5P_TKT', ['tktC2']), ('F6P_GAP_TAL', ['talC3']), ('S7P_E4P_TAL', ['tal']), ('FBP', ['FBP_Vmax']),
                      ('PPC', ['PPC_Vmax']), ('PCK', ['PCK_Vmax']), ('PPS', ['PPS_Vmax']), ('MAD', ['MAD_Vmax']),
                      ('PDH', ['PDH_Vmax']), ('GLT', ['GLT_Vmax']), ('ACN_1', ['ACN_1_Vmax']), ('ACN_2', ['ACN_2_Vmax']),
                      ('ICD', ['icd']), ('LPD', ['LPD_Vmax']), ('SK', ['SK_Vmax']), ('SDH', ['SDH_Vmax']), ('FUMA', ['FUMA_Vmax']),
                      ('MQO', ['MQO_Vmax']), ('MDH', ['MDH_Vmax']), ('ACEA', ['ACEA_Vmax']), ('ACEB', ['ACEB_Vmax']),
                      ('EDD', ['EDD_Vmax']), ('EDA', ['EDA_Vmax']), ('NADH_req', ['NADH_req_Vmax']), ('ATP_syn', ['ATP_syn_Vmax']),
                      ('ACK', ['ACK_Vmax']), ('ACS', ['ACS_Vmax']), ('PTA', ['PTA_Vmax']), ('MYTBO', ['MYTBO_Vmax']),
                      ('SQR', ['SQR_Vmax']), ('NDHII', ['NDHII_Vmax']), ('GROWTH', ['GROWTH_Vmax']), ('ATP_MAINTENANCE', ['ATP_MAINTENANCE_Vmax']),
                      ('XCH_GLC', ['XCH_GLC_Vmax']), ('PIT', ['PIT_Vmax']), ('XCH_P', ['XCH_P_Vmax']), ('XCH_ACE1', ['XCH_ACE1_Vmax']),
                      ('XCH_ACE2', ['XCH_ACE2_Vmax'])])

    model = load_kinetic_model(sbmlFile, mapParamReacs)

    objFunc = build_evaluation_function("targetFlux", ["_ACE_OUT"])
    simulProblem= KineticSimulationProblem(model, tSteps=[0, KineticConfigurations.STEADY_STATE_TIME])
    result = kinetic_strain_optim(simulProblem, objFunc=objFunc, levels=LEVELS, criticalParameters=['ATP_MAINTENANCE_Vmax', 'GROWTH_Vmax', 'NDHII_Vmax', 'PIT_Vmax', 'eiicbP', 'ei', 'eiP', 'eiia'], isMultiProc=isMultiProc, resultFile=fileRes, initPopFile=fileLastGen)
    result.print()


def ko_CoreColi(isMultiProc=False, size=1):
    SBML_FILE = '../../../examples/models/CoreEcoli_fixed.xml'
    basePath = "C:/Users/BISBII/Results/"
    fileRes = basePath + "optim222_CoreColi_succ_ko_sizex" + str(size) + ".csv"
    map = {}
    for i in range(1,92):
        reaction = 'V'+str(i)
        param = 'vMax_R' + str(i)
        map[reaction] = [param]

    model = load_kinetic_model(SBML_FILE, map)
    objFunc = build_evaluation_function("targetFlux", ["V57"])
    simulProblem = KineticSimulationProblem(model, tSteps=[0, KineticConfigurations.STEADY_STATE_TIME])
    simulProblem.simulate()

    #kinetic_strain_optim(simulProblem, objFunc=objFunc, type=optimType.REACTION_KO, isMultiProc=isMultiProc, candidateSize=size, resultFile=fileRes)


    '''
    over = {'vMax_R11': 0}
    ko = OrderedDict(over)
    from optimModels.simulation.override_simul_problem import OverrideKineticSimulProblem
    override = OverrideKineticSimulProblem(ko)
    res = simulProblem.simulate(overrideSimulProblem=override)
    res.print()
    '''

mapa={"R_PDH":["rmaxPDH"],"R_TALA": ["rmaxTALA","KTALAeq"], "R_PGMT": ["rmaxPGMT"],"R_TKT1": ["rmaxTKT1"],"R_TKT2": ["rmaxTKT2"],"R_PTS":["rmaxPTS"],"R_PGI":["rmaxPGI"],"R_PYK":["rmaxPYK"]}

def ko_model(isMultiProc=False, size=4):
    sbmlFile = 'C:/Users/BiSBII/Downloads/chassagnole2002_mMs.xml'
    fileRes ="Optim/optim_chassa_pdh_ko_"+str(size)+".csv"


    model = load_kinetic_model(sbmlFile,mapa)
    print(len(model.reactions))

    #objFunc = build_evaluation_function("targetFlux", ["R_TKT1"])
    #simulProblem = KineticSimulationProblem(model, tSteps=[0, KineticConfigurations.STEADY_STATE_TIME])

    #sim = simulProblem.simulate()
    #sim.print()
    #results = kinetic_strain_optim(simulProblem, objFunc=objFunc, levels=None, type= optimType.REACTION_KO, isMultiProc=isMultiProc, candidateSize=size, resultFile=fileRes)
    #for result in results:
    #    result.print()

if __name__ == '__main__':
    import time
    size = 1
    #t1 = time.time()
    start_time = time.time()
    #ko_chassagnole(isMultiProc=True, size=size)
    #ko_CoreColi(isMultiProc=False, size=size)
    ko_model()
    print("--- %s seconds ---" % (time.time() - start_time), 'end')
    #ko_millard(size=7)
    #ko_jahan()
