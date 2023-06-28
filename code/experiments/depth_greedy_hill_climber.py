from ..algorithms.greedy import Greedy
from ..algorithms.hill_climber import Hill_climber
from ..algorithms.depth_first import DepthFirst
from ..classes.protein import Protein
from ..algorithms.simulated_annealing import SimulatedAnnealing

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path = r"C:\Users\alaya\Documents\Data_experiments"

# create dataframe to store results
df_exp_depth_greedy_hill = pd.DataFrame()
df_exp_depth_greedy_hill = pd.DataFrame(columns = ['Depth_first_score', 'Depth_first_states', 'Greedy_score', 'Greedy_states', 'Hill_climber_scores', 'Hill_climber_states', 'Sim_anneal_scores', 'Sim_anneal_states'])

# initialize lists and variables
greedy_scores = []
hill_climber_scores = []
sim_anneal_scores = []

greedy_states = 0
hill_climber_states = 0
depth_first_states = 531441
sim_anneal_states = 0

# run depth first with test protein
test_protein = Protein("HHPHHHPHPHHHPH")

depth_first = DepthFirst(test_protein, 2)
depth_first.run(P_pruning=False, directions_pruning=False)
print('Depth First completed')

df_exp_depth_greedy_hill.to_csv(path_or_buf=fr"{path}\df_exp_depth_greedy_hill_complete")

# make a dictionary with the splits as keys and the number of iterations per 'before' as the value
dict_splits = {1: 10, 2: 5, 3: 4, 4: 3}
total_states = 0
# run greedy with depth with different split sizes and before sizes with test protein
for key, value in dict_splits.items():
    for before in range(key):
        for i in range(value):
            test_protein = Protein("HHPHHHPHPHHHPH")

            greedy_protein = Greedy(test_protein, 2, splits=key, before=before)
            greedy_protein.run()

            greedy_scores.append(greedy_protein.protein.score)
            greedy_states += greedy_protein.states

df_exp_depth_greedy_hill.to_csv(path_or_buf=fr"{path}\df_exp_depth_greedy_hill_complete")
# run hill climber with test protein
n = 50

for i in range(n):
    test_protein = Protein("HHPHHHPHPHHHPH")

    hill_climber_protein = Hill_climber(test_protein, 2)
    hill_climber_protein.run_i_iterations(test_protein, iterations=500, bonds=1)

    hill_climber_scores.append(hill_climber_protein.protein.score)
    hill_climber_states += hill_climber_protein.states

print('Hillclimber completed')

df_exp_depth_greedy_hill.to_csv(path_or_buf=fr"{path}\df_exp_depth_greedy_hill_complete")
# run simulated annealing with test protein
n = 50
for i in range(n):
    test_protein = Protein("HHPHHHPHPHHHPH")

    simanneal_protein = SimulatedAnnealing(test_protein, start_n=10, folded=False, dimensions = 2, temperature = 10)
    simanneal_protein.run_i_iterations(test_protein, iterations = 500, bonds=10)

    sim_anneal_scores.append(simanneal_protein.protein.score)
    sim_anneal_states += simanneal_protein.states

print('Simulated Annealing completed')

df_exp_depth_greedy_hill.to_csv(path_or_buf=fr"{path}\df_exp_depth_greedy_hill_complete")

# store all results in dataframe and save dataframe on computer
df_exp_depth_greedy_hill['Depth_first_score'] = [depth_first.protein.score]
df_exp_depth_greedy_hill['Depth_first_states'] = [depth_first_states]
                                                      
df_exp_depth_greedy_hill['Greedy_score'] = [min(greedy_scores)]
df_exp_depth_greedy_hill['Greedy_states'] = [greedy_states]

df_exp_depth_greedy_hill['Hill_climber_scores'] = [min(hill_climber_scores)]
df_exp_depth_greedy_hill['Hill_climber_states'] = [hill_climber_states]

df_exp_depth_greedy_hill['Sim_anneal_scores'] = [min(sim_anneal_scores)]
df_exp_depth_greedy_hill['Sim_anneal_states'] = [sim_anneal_states]

df_exp_depth_greedy_hill.to_csv(path_or_buf=fr"{path}\df_exp_depth_greedy_hill_complete")

# -------------------- VISUALIZE --------------------
df_exp_depth_greedy_hill = pd.read_csv(fr"{path}\df_exp_depth_greedy_hill_complete")

# create combined barchart of states and scores for all algorithms
states = [int(df_exp_depth_greedy_hill['Depth_first_states']), int(df_exp_depth_greedy_hill['Greedy_states']), int(df_exp_depth_greedy_hill['Hill_climber_states']), int(df_exp_depth_greedy_hill['Sim_anneal_states'])]
scores = [int(df_exp_depth_greedy_hill['Depth_first_score']), int(df_exp_depth_greedy_hill['Greedy_score']), int(df_exp_depth_greedy_hill['Hill_climber_scores']), int(df_exp_depth_greedy_hill['Sim_anneal_scores'])]

fig, ax1 = plt.subplots()

x_array = np.arange(4)
width = float(0.40)

for i in range(len(x_array)):
    ax1.bar(x_array[i] - 0.2, states[i], width=width, color='orange')
    ax1.set_xlabel("Algorithms")
    ax1.set_ylabel("states")

    ax2 = ax1.twinx()
    ax2.bar(x_array[i] + 0.2, scores[i], width=width, color='blue')
    ax2.set_ylabel("Scores")
    ax2.set_ylim(min(scores), 0)
    ax2.invert_yaxis()

ax1.legend(["States"], loc=2)
ax2.legend(["Scores"], loc=1)
plt.xticks(x_array, ['Depth First', 'Greedy with depth', 'Hill Climber', 'Simulated Annealing'])
plt.title('Comparison of scores and states of short protein in 2D')

plt.savefig(fr"{path}\Depth_First_vs_Greedy_vs_Hill_Climber_vs_Simulated_Annealing_complete.pdf")
plt.show()