import sys
import pandas as pd
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(suppress=True)
import matplotlib.pyplot as plt
pd.options.display.width = 0
pd.set_option('display.max_rows', None)

def sol_str_to_list (str_sol):
    str_sol = str_sol.rstrip()
    str_sol = str_sol[1:-1].split(sep = ', ')
    str_sol = [x[1:-1] for x in str_sol]
    list_sol = str_sol
    return list_sol

def save_to_csv(csv_address, line, header = False):
    file = open(csv_address, "a+")
    if header:
        file.write(line)
    else:
        file.write("\n"+line)
    file.close()

# Upload data
results_dir = "C:Users/theze/Documentos/Work/CEB_Fellowship/Work/Results/Succ_Optim/Yeast_Succ/"
csv_name = "wyield1500_AERO_GECKO_KO_SUCC_max7.csv"
data = pd.read_csv(csv_name, sep = ";", skiprows = 2, skipfooter = 103, engine = "python")


# # Evolution plot
# data_best = data.sort_values('Fitness', ascending = False).drop_duplicates(['Generation']).sort_index()
# scplot = data_best.plot.scatter(x = "Generation", y = "Fitness")
# plt.title('{}'.format(csv_name))
# plt.show()

# Find most used proteins
data = data.sort_values('Fitness', ascending = False).drop_duplicates(['Reactions']).sort_index()
data = data[data['Fitness'] > 1].Reactions

sol_counting = {}
for solu in data:
    for prot in sol_str_to_list(solu):
        if prot in sol_counting:
            sol_counting[prot] += 1
        else:
            sol_counting[prot] = 1

sol_counting = [(k, v) for k, v in sol_counting.items()]
sol_counting = sorted(sol_counting, key = lambda x: x[1], reverse = True)

from Bio import SeqIO
import urllib
for protein in sol_counting:
    handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/{}.xml".format(protein[0]))
    record = SeqIO.read(handle, "uniprot-xml")
    name = record.description
    print("{}: {} - {}".format(protein[1], protein[0], name))

cytochromes = ["P00128", "P07256", "P08067", "P07257", "P37299", "P00127", "P07143", "P00163", "P08525", "P22289", "P37298", "P33421", "P32799 ", "P53239", "P00045", "P00175", "Q01519", "P04037", "P19516", "Q3E731", "Q04935", "P54781", "P38824", "P47081"]
glycine_cleavage_complex = ["P49095", "P48015", "P39726"]
proteins_to_exclude = cytochromes + glycine_cleavage_complex

# # Solution cluster data
# good_sols = data.sort_values('Fitness', ascending = False).drop_duplicates(['Reactions'])
# good_sols_nozero = good_sols[good_sols['Fitness'] >= 1]
# gsnz_reac_only = good_sols_nozero.Reactions
#
# all_prots = []
# for sol in gsnz_reac_only:
#     sol = sol_str_to_list(sol)
#     all_prots += sol
# all_prots = sorted(list(set(all_prots)))
# max_value = len(all_prots)
#
# values = np.zeros((max_value, max_value))
#
# for index, row in good_sols_nozero.iterrows():
#     sol = sol_str_to_list(row['Reactions'])
#     for indx in range(len(sol)):
#         for indx2 in range(indx + 1, len(sol)):
#             value1 = all_prots.index(sol[indx])
#             value2 = all_prots.index(sol[indx2])
#             fit_value = (float("{0:.2f}".format(row['Fitness'])), 0.00)[float("{0:.2f}".format(row['Fitness'])) < 0.01]
#             values[value1][value2] += fit_value
#             values[value2][value1] += fit_value
#
#
# filename = 'C:/Users/theze/PycharmProjects/optimmodels_jms/src/optimModels/unittests/gephi_files/{}'.format(csv_name)
# header = ';' + str(all_prots)[1:-1].replace(',', ';').replace(" ", "").replace("'", "")
# save_to_csv(filename, header, header = True)
#
# vls_indx = range(max_value)
# for val_indx in vls_indx:
#     line_0 = str(all_prots[val_indx]) + ';' + " ".join(str(values[val_indx])[1:-1].replace("\n", "").replace("0. ", "0.00 ").rstrip().lstrip().split()).replace(" ", ";")
#     save_to_csv(filename, line_0)