from ..algorithms.greedy import Greedy
from ..algorithms.hill_climber import Hill_climber
from ..algorithms.depth_first import DepthFirst
from ..classes.protein import Protein
from ..algorithms.simulated_annealing import SimulatedAnnealing

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path = r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data"

# create dataframe to store results
df_exp_depth_greedy_hill = pd.DataFrame()
df_exp_depth_greedy_hill = pd.DataFrame(columns = ['Depth_first_score', 'Depth_first_iterations', 'Greedy_score', 'Greedy_iterations', 'Hill_climber_scores', 'Hill_climber_iterations', 'Sim_anneal_scores', 'Sim_anneal_iterations'])

# initialize lists and variables
greedy_scores = []
hill_climber_scores = []
sim_anneal_scores = []

greedy_iterations = 0
hill_climber_iterations = 0
depth_first_iterations = 1
sim_anneal_iterations = 0

# run depth first with test protein
test_protein = Protein("HHPHHHPHPHHHPH")

depth_first = DepthFirst(test_protein, 2)
depth_first.run(P_pruning=True, directions_pruning=True)
print('Depth First completed')

# make a dictionary with the splits as keys and the number of iterations per 'before' as the value
dict_splits = {1: 2000, 2: 1000, 3: 667, 4: 500, 5: 400}

# run greedy with depth with different split sizes and before sizes with test protein
for key, value in dict_splits.items():
    for before in range(key):
        for i in range(value):
            test_protein = Protein("HHPHHHPHPHHHPH")

            greedy_protein = Greedy(test_protein, 2, splits=key, before=before)
            greedy_protein.run()

            greedy_scores.append(greedy_protein.protein.score)
            greedy_iterations += 1

# run hill climber with test protein
n = 50

for i in range(n):
    test_protein = Protein("HHPHHHPHPHHHPH")

    hill_climber_protein = Hill_climber(test_protein, 2)
    hill_climber_protein.run_i_iterations(test_protein, iterations=500, bonds=1)

    hill_climber_scores.append(hill_climber_protein.protein.score)
    hill_climber_iterations += 1

print('Hillclimber completed')

# run simulated annealing with test protein
for i in range(n):
    for j in range(1, 11):
        test_protein = Protein("HHPHHHPHPHHHPH")

        # CHANGE TEMP BASED ON IDEAL SIMULATED ANNEALING
        simanneal_protein = SimulatedAnnealing(test_protein, start_n=j, folded=False, dimensions=2, temperature=20)
        simanneal_protein.run_i_iterations(test_protein, iterations=500, bonds=j)

        sim_anneal_scores.append(simanneal_protein.protein.score)
        sim_anneal_iterations += 1

print('Simulated Annealing completed')

# store all results in dataframe and save dataframe on computer
df_exp_depth_greedy_hill['Depth_first_score'] = [depth_first_scores_list]
df_exp_depth_greedy_hill['Depth_first_iterations'] = [depth_first_iterations_list]
                                                      
df_exp_depth_greedy_hill['Greedy_score'] = [greedy_scores_list]
df_exp_depth_greedy_hill['Greedy_iterations'] = [greedy_iterations_list]

df_exp_depth_greedy_hill['Hill_climber_scores'] = [hill_climber_scores_list]
df_exp_depth_greedy_hill['Hill_climber_iterations'] = [hill_climber_iterations_list]

df_exp_depth_greedy_hill['Sim_anneal_scores'] = [sim_anneal_scores_list]
df_exp_depth_greedy_hill['Sim_anneal_iterations'] = [sim_anneal_iterations_list]

df_exp_depth_greedy_hill.to_csv(path_or_buf=fr"{path}\df_exp_depth_greedy_hill_complete")

# -------------------- VISUALIZE --------------------
df_exp_depth_greedy_hill = pd.read_csv(fr"{path}\df_exp_depth_greedy_hill_complete")

# create combined barchart of iterations and scores for all algorithms
iterations = [int(df_exp_depth_greedy_hill['Depth_first_iterations']), int(df_exp_depth_greedy_hill['Greedy_iterations']), int(df_exp_depth_greedy_hill['Hill_climber_iterations']), int(df_exp_depth_greedy_hill['Sim_anneal_iterations'])]
scores = [int(df_exp_depth_greedy_hill['Depth_first_score']), int(df_exp_depth_greedy_hill['Greedy_score']), int(df_exp_depth_greedy_hill['Hill_climber_scores']), int(df_exp_depth_greedy_hill['Sim_anneal_scores'])]

fig, ax1 = plt.subplots()

x_array = np.arange(4)
width = float(0.40)

for i in range(len(x_array)):
    ax1.bar(x_array[i] - 0.2, iterations[i], width=width, color='orange')
    ax1.set_xlabel("Algorithms")
    ax1.set_ylabel("Iterations")

    ax2 = ax1.twinx()
    ax2.bar(x_array[i] + 0.2, scores[i], width=width, color='blue')
    ax2.set_ylabel("Scores")
    ax2.set_ylim(min(scores), 0)
    ax2.invert_yaxis()

ax1.legend(["Iterations"], loc=2)
ax2.legend(["Scores"], loc=1)
plt.xticks(x_array, ['Depth First', 'Greedy with depth', 'Hill Climber', 'Simulated Annealing'])
plt.title('Comparison of scores and iterations of short protein in 2D')

plt.savefig(fr"{path}\Depth_First_vs_Greedy_vs_Hill_Climber_vs_Simulated_Annealing_complete.pdf")
plt.show()