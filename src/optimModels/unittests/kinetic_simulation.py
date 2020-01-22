from optimModels.model.kineticModel import load_kinetic_model
from optimModels.simulation.run import kinetic_simulation
import pandas as pd
from random import randint

class Params :
    def __init__(self):
        file = '../../../examples/data/alternativeParamSets_toyModel.csv'
        aux= pd.read_csv(file)
        self.params = aux.iloc[:,1:]
        self.params.index = list(aux.iloc[:, 0])

    def get_params_set(self):
        selected = randint(0,self.params.shape[1]-1)

        data = self.params.iloc[:,selected]
        print (data.to_dict)
        return data.to_dict()

# Optimization:
# model: TinyModel_RHS (EPFL model)
# Objective: Validate in the toy model works with several combinations of parameters
def toy_model():
    SBML_FILE = '../../../examples/models/TinyModel_RHS.xml'

    mapParamReacs = "(vMax\w*)" # regular expression to identify the vmax parameter

    model = load_kinetic_model(SBML_FILE, [])

    file = '../../../examples/data/alternativeParamSets_toyModel.csv'
    aux = pd.read_csv(file)
    params = aux.iloc[:, 1:]
    params.index = list(aux.iloc[:, 0])

    for i in range(1000):
        data = params.iloc[:, i]
        res = kinetic_simulation(model, parameters = data.to_dict(), factors = None, time = 1e9)

        print( str(i)+ " -->" + str(res.solverStatus))

 #   for  k,v in model.get_parameters().items():
 #       print(k + " --> " + str(v))
 #   for r, f in res.get_fluxes_distribution().items():
 #       print (r + " --> " +str(f))

  #  for m, c in res.get_steady_state_concentrations().items():
  #      print(m + " --> " + str(c))


def toy_model2():
    SBML_FILE = '../../../examples/models/C_Mutant_0.xml'

    model = load_kinetic_model(SBML_FILE, [])

    res = kinetic_simulation(model, time=1e9)
    res.print()
    #print (res.get_fluxes_distribution())


# Optimization:
# model: Jahan2016_chemostat_fixed
# Objective: Test the under/over expression of some vMax parameters
def kinetic_sim():
    SBML_FILE = '../../../examples/models/Jahan2016_chemostat_fixed.xml'
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

    model = load_kinetic_model(SBML_FILE, mapParamReacs)
    factors = None
    factors = {'CS': 32, 'v6PGDH_max': 0.0625, 'vEdd_max': 8, 'ICDH': 0.125, 'MS': 8, 'vPta_max': 0, 'SDH': 0.03125, 'vPTS4_max': 32}
    factors = {'ICDH': 0.03125, 'vAck_max': 0, 'SDH': 0.0625, 'vPTS4_max': 32}
    res = kinetic_simulation(model, parameters = None, factors = factors, time = 1e9)

    for r, f in res.get_fluxes_distribution().items():
        print (r + " --> " +str(f))

    for m, c in res.get_steady_state_concentrations().items():
        print(m + " --> " + str(c))

def CoreEcoli():
    SBML_FILE = '../../../examples/models/CoreEcoli_fixed.xml'
    map = {}
    for i in range(1,92):
        reaction = 'V'+str(i)
        param = 'vMax_R' + str(i)
        map[reaction] = [param]

    print(map)

    model = load_kinetic_model(SBML_FILE, map)
    #params = list(model.get_parameters().keys())[1:-2]

    #alternative_parameters = {}

    #fic = open("../../../examples/models/CoreEcoli_AlternativeParameters.csv", 'r')
    #lines = fic.readlines()
    #for i in range(len(params)):
    #    alternative_parameters[params[i]] = lines[i].strip().split(',')




    #from framed import time_course
    #T, X = time_course(model, time=100)
    #print(X)

    from optimModels.simulation.simul_problems import KineticSimulationProblem
    simulProb = KineticSimulationProblem(model, timeout=None)
    #print(simulProb.get_time_steps())
    res = simulProb.simulate()
    res.print()

def get_reactions_names():
    model_file = '../../../examples/models/iJO1366.xml'
    names = '../../../examples/models/reactionIDs.csv'

    model = load_kinetic_model(model_file, {})

    import pandas as pd
    df = pd.read_csv(names, header=0)
    ids = df['ReactionName'].tolist()
    dic = {}
    for r in model.reactions.values():
        dic[r.id[2:]] = r.name
    reaction_names = []
    for r in ids:
        if r in dic.keys():
            name = dic[r]
        else: name = '---'
        reaction_names.append(name)
    new_df = pd.DataFrame({'ids':ids, 'names':reaction_names})
    new_df.to_csv('reaction_ids_and_names.csv')
    #print(reaction_names)



if __name__ == '__main__':
    CoreEcoli()
    #get_reactions_names()
    #toy_model2()
