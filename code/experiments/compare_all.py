from code.classes.protein import Protein

from code.algorithms.randomise import *
from code.algorithms.greedy import Greedy
from code.algorithms.depth_first import DepthFirst
from code.algorithms.important_parts import ImportantParts
from code.algorithms.hill_climber import Hill_climber
from code.algorithms.simulated_annealing import SimulatedAnnealing

import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path = r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data"

# create dataframe to store results
df_short_scores = pd.DataFrame()
df_long_scores = pd.DataFrame()

test_protein_short = Protein("HHPHHHPHPHHHPH")
test_protein_long = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

# ---------------------------- Random ----------------------------
# run random algorithm for short and long protein and store results
start = time.time()

random_protein_short, score_list_short, best_score_short = random_reassignment(test_protein_short, 3, k=1000000)

end = time.time()
random_time_short = end - start

start = time.time()

random_protein_long, score_list_long, best_score_long = random_reassignment(test_protein_long, 3, k=1000000)

end = time.time()
random_time_long = end - start

df_short_scores["Random score"] = [random_protein_short.score]
df_short_scores["Random time"] = [random_time_short]

df_long_scores["Random score"] = [random_protein_long.score]
df_long_scores["Random time"] = [random_time_long]

print('Random done!')

df_short_scores.to_csv(path_or_buf=fr"{path}\df_short_scores_complete")
df_long_scores.to_csv(path_or_buf=fr"{path}\df_long_scores_complete")
# ---------------------------- Greedy ----------------------------
# run greedy algorithm for short and long protein and store results
greedy_scores_short = []
greedy_scores_long = []

start = time.time()

n = 25000

for i in range(n):
    test_protein_short = Protein("HHPHHHPHPHHHPH")

    greedy_protein_short = Greedy(test_protein_short, 3, splits=1)
    greedy_protein_short.run()
    greedy_scores_short.append(greedy_protein_short.protein.score)

end = time.time()
greedy_time_short = end - start

start = time.time()
for i in range(n):
    test_protein_long = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

    greedy_protein_long = Greedy(test_protein_long, 3, splits=1)
    greedy_protein_long.run()
    greedy_scores_long.append(greedy_protein_long.protein.score)

end = time.time()
greedy_time_long = end - start

df_short_scores["Greedy score"] = [min(greedy_scores_short)]
df_short_scores["Greedy time"] = [greedy_time_short]

df_long_scores["Greedy score"] = [min(greedy_scores_long)]
df_long_scores["Greedy time"] = [greedy_time_long]

print('Greedy done!')

df_short_scores.to_csv(path_or_buf=fr"{path}\df_short_scores_complete")
df_long_scores.to_csv(path_or_buf=fr"{path}\df_long_scores_complete")
# ---------------------------- Greedy with depth search ----------------------------
# run greedy with depth algorithm for short and long protein and store results
greedy_beam_scores_short = []
greedy_beam_scores_long = []

start = time.time()

n = 500

for i in range(n):
    for before in range(4):
        test_protein_short = Protein("HHPHHHPHPHHHPH")

        greedy_beam_protein_short = Greedy(test_protein_short, 3, splits=4, before=before)
        greedy_beam_protein_short.run()
        greedy_beam_scores_short.append(greedy_beam_protein_short.protein.score)

end = time.time()
greedy_beam_time_short = end - start

start = time.time()

for i in range(n):
    for before in range(4):
        test_protein_long = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

        greedy_beam_protein_long = Greedy(test_protein_long, 3, splits=4, before=before)
        greedy_beam_protein_long.run()
        greedy_beam_scores_long.append(greedy_beam_protein_long.protein.score)

end = time.time()
greedy_beam_time_long = end - start

df_short_scores["Greedy score with depth"] = [min(greedy_beam_scores_short)]
df_short_scores["Greedy depth time"] = [greedy_beam_time_short]

df_long_scores["Greedy score with depth"] = [min(greedy_beam_scores_long)]
df_long_scores["Greedy depth time"] = [greedy_beam_time_long]

print('Greedy with depth done!')

df_short_scores.to_csv(path_or_buf=fr"{path}\df_short_scores_complete")
df_long_scores.to_csv(path_or_buf=fr"{path}\df_long_scores_complete")

# ---------------------------- Depth First with P and directions pruning ----------------------------
# run depth first algorithm with P pruning and directions pruning for short protein and store results
start = time.time()

depth_first_protein = DepthFirst(test_protein_short, 3)
depth_first_protein.run(P_pruning=True, directions_pruning=True)

end = time.time()
depth_first_time = end - start

df_short_scores["Depth First score"] = [depth_first_protein.protein.score]
df_short_scores["Depth First time"] = [depth_first_time]

print('Depth First done!')

df_short_scores.to_csv(path_or_buf=fr"{path}\df_short_scores_complete")
df_long_scores.to_csv(path_or_buf=fr"{path}\df_long_scores_complete")
# ---------------------------- Depth First - Important Parts ----------------------------
# run important parts algorithm for short and long protein and store results
important_parts_scores_short_on_P = []
important_parts_scores_long_on_P = []

start = time.time()
important_parts_protein_short_on_P = ImportantParts(test_protein_short, 3)
important_parts_protein_short_on_P.run(iterations=1000, split_on_P=True, split_on_size=False)

end = time.time()
important_parts_time_short_on_P = end - start

start = time.time()
important_parts_protein_long_on_P = ImportantParts(test_protein_long, 3)
important_parts_protein_long_on_P.run(iterations=1000, split_on_P=True, split_on_size=False)

end = time.time()
important_parts_time_long_on_P = end - start

df_short_scores["Important Parts score"] = [important_parts_protein_short_on_P.protein.score]
df_short_scores["Important Parts time"] = [important_parts_time_short_on_P]

df_long_scores["Important Parts score"] = [important_parts_protein_long_on_P.protein.score]
df_long_scores["Important Parts time"] = [important_parts_time_long_on_P]

print('Important Parts done!')

df_short_scores.to_csv(path_or_buf=fr"{path}\df_short_scores_complete")
df_long_scores.to_csv(path_or_buf=fr"{path}\df_long_scores_complete")

# ---------------------------- Hill Climber ----------------------------
# run hill climber algorithm for short and long protein and store results
start = time.time()

test_protein_short = Protein("HHHHHHHH")
hill_climber_protein_short = Hill_climber(test_protein_short, dimensions = 3, folded = False)
hill_climber_protein_short.experiment(test_protein_short, iterations = 500, sample_size = 50, max_n = 1, plot = False)

end = time.time()
hill_climber_time_short = end - start

start = time.time()
test_protein_long = Protein("HHHHHHHH")
hill_climber_protein_long = Hill_climber(test_protein_long, dimensions = 3, folded = False)
hill_climber_protein_long.experiment(test_protein_long, iterations = 500, sample_size = 50, max_n = 1, plot = False)

end = time.time()
hill_climber_time_long = end - start

df_short_scores["Hill climber score"] = [hill_climber_protein_short.protein.score]
df_short_scores["Hill climber time"] = [hill_climber_time_short]

df_long_scores["Hill climber score"] = [hill_climber_protein_long.protein.score]
df_long_scores["Hill climber time"] = [hill_climber_time_long]

print('Hill Climber done!')

df_short_scores.to_csv(path_or_buf=fr"{path}\df_short_scores_complete")
df_long_scores.to_csv(path_or_buf=fr"{path}\df_long_scores_complete")

# ---------------------------- Simulated Annealing ----------------------------
# run simulated annealing algorithm for short and long protein and store results
# fill in the right parameters
start = time.time()

test_protein_short = Protein("HHHHHHHHHH")
simanneal_protein_short = SimulatedAnnealing(test_protein_short, start_n = 10, folded = False, dimensions = 3,  temperature = 20)
simanneal_protein_short.run_i_iterations(test_protein_short, iterations=1000, bonds=10, sim_annealing=True)

end = time.time()
sim_anneal_short_time = end - start

start = time.time()

test_protein_long = Protein("HHHHHHHHHH")
simanneal_protein_long = SimulatedAnnealing(test_protein_long, start_n = 10, folded = False, dimensions = 3,  temperature = 20)
simanneal_protein_long.run_i_iterations(test_protein_long, iterations=1000, bonds=10, sim_annealing=True)

end = time.time()
sim_anneal_long_time = end - start

df_short_scores["Simulated Annealing score"] = [simanneal_protein_short.protein.score]
df_short_scores["Simulated Annealing time"] = [sim_anneal_short_time]

df_long_scores["Simulated Annealing score"] = [simanneal_protein_long.protein.score]
df_long_scores["Simulated Annealing time"] = [sim_anneal_long_time]

print('Simulated Annealing done!')

df_short_scores.to_csv(path_or_buf=fr"{path}\df_short_scores_complete")
df_long_scores.to_csv(path_or_buf=fr"{path}\df_long_scores_complete")

# -------------------- VISUALIZE SHORT COMPARISON--------------------
df_short_scores = pd.read_csv(fr"{path}\df_short_scores_complete")

# create a combined barchart of the scores and the time for all algorithms for the short protein
time_short = [int(df_short_scores ['Random time'].iloc[0]), int(df_short_scores['Greedy time'].iloc[0]), int(df_short_scores['Greedy depth time'].iloc[0]), int(df_short_scores['Depth First time'].iloc[0]), int(df_short_scores['Important Parts time'].iloc[0]), int(df_short_scores['Hill climber time'].iloc[0]), int(df_short_scores['Simulated Annealing time'].iloc[0])]
scores_short = [int(df_short_scores['Random score'].iloc[0]), int(df_short_scores['Greedy score'].iloc[0]), int(df_short_scores['Greedy score with depth'].iloc[0]), int(df_short_scores['Depth First score'].iloc[0]), int(df_short_scores['Important Parts score'].iloc[0]), int(df_short_scores['Hill climber score'].iloc[0]), int(df_short_scores['Simulated Annealing score'].iloc[0])]

fig, ax1 = plt.subplots()

x_array = np.arange(7)
width = float(0.40)

for i in range(len(x_array)):
    ax1.bar(x_array[i] - 0.2, time_short[i], width=width, color='orange')
    ax1.set_xlabel("Algorithms")
    ax1.set_ylabel("Time")
    
    ax2 = ax1.twinx()

    ax2.bar(x_array[i] + 0.2, scores_short[i], width=width, color='blue')
    ax2.set_ylabel("Scores")
    ax2.set_ylim(min(scores_short), 0)
    ax2.invert_yaxis()

plt.xticks(x_array, ['Random', 'Greedy', 'Greedy with depth', 'Depth First', 'Important Parts', 'Hill Climber', 'Simulated Annealing'])
ax1.legend(["Time (s)"], loc=2)
ax2.legend(["Scores"], loc=1)

plt.title('Comparison of all algorithms for short protein')
plt.savefig(fr"{path}\Full_comparison_short_complete.pdf")
plt.show()

# -------------------- VISUALIZE LONG COMPARISON--------------------
df_long_scores = pd.read_csv(fr"{path}\df_long_scores_complete")

# create a combined barchart of the scores and the time for all algorithms for the long protein
time_long = [int(df_long_scores['Random time'].iloc[0]), int(df_long_scores['Greedy time'].iloc[0]), int(df_long_scores['Greedy depth time'].iloc[0]), int(df_long_scores['Important Parts time'].iloc[0]), int(df_long_scores['Hill climber time'].iloc[0]), int(df_long_scores['Simulated Annealing time'].iloc[0])]
scores_long = [int(df_long_scores['Random score'].iloc[0]), int(df_long_scores['Greedy score'].iloc[0]), int(df_long_scores['Greedy score with depth'].iloc[0]), int(df_long_scores['Important Parts score'].iloc[0]), int(df_long_scores['Hill climber score'].iloc[0]), int(df_long_scores['Simulated Annealing score'].iloc[0])]

fig, ax1 = plt.subplots()

x_array = np.arange(6)
width = float(0.40)

for i in range(len(x_array)):
    ax1.bar(x_array[i] - 0.2, time_long[i], width=width, color='orange')
    ax1.set_xlabel("Algorithms")
    ax1.set_ylabel("Time")

    ax2 = ax1.twinx()

    ax2.bar(x_array[i] + 0.2, scores_long[i], width=width, color='blue')
    ax2.set_ylabel("Scores")
    ax2.set_ylim(min(scores_long), 0)
    ax2.invert_yaxis()

plt.xticks(x_array, ['Random', 'Greedy', 'Greedy with depth', 'Important Parts', 'Hill Climber', 'Simulated Annealing'])
ax1.legend(["Time (s)"], loc=2)
ax2.legend(["Scores"], loc=1)

plt.title('Comparison of all algorithms for long protein')
plt.savefig(fr"{path}\Full_comparison_long_complete.pdf")
plt.show()