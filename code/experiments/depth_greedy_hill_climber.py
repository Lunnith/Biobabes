from ..algorithms.greedy import Greedy
from ..algorithms.hill_climber import Hill_climber
from ..algorithms.depth_first import DepthFirst

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_exp_depth_greedy_hill = pd.DataFrame()

greedy_scores = []
hill_climb_scores = []
sim_anneal_scores = []

greedy_iterations = 0
hill_climber_iterations = 0
depth_first_iterations = 1
sim_anneal_iterations = 0

test_protein = Protein("HHPHHHPHPHHHPH")

depth_first = DepthFirst(test_protein, 2)
depth_first.run(P_pruning=True, directions_pruning=True)

# CHANGE
n = 1

for i in range(n):
    for i in range(1, 5):
        test_protein = Protein("HHPHHHPHPHHHPH")

        before_list = []

        for number in range(i):
            if number < i:
                before_list.append(number)

        for before in before_list:
            greedy_protein = Greedy(test_protein, 2, splits=i, before=before)
            greedy_protein.run()

            greedy_scores.append(greedy_protein.protein.score)
            greedy_iterations += 1
    
for i in range(n):
    test_protein = Protein("HHPHHHPHPHHHPH")

    hill_climber_protein = Hill_climber(test_protein, 2)
    hill_climber_protein.run_i_iterations(test_protein, iterations=500, bonds=1)

    hill_climber_scores.append(hill_climber_protein.protein.score)
    hill_climber_iterations += 1

for i in range(n):
    for i in range(1, 10):
        test_protein = Protein("HHPHHHPHPHHHPH")

        # CHANGE TEMP BASED ON IDEAL SIMULATED ANNEALING
        simanneal_protein = SimulatedAnnealing(test_protein, i, 2, temperature=20)
        simanneal_protein.run_i_iterations(test_protein, iterations=500, bonds=i)

        sim_anneal_scores.append(simanneal_protein.protein.score)
        sim_anneal_iterations += 1

df_exp_depth_greedy_hill = pd.DataFrame()

df_exp_depth_greedy_hill['depth_first_score'] = depth_first.protein.score
df_exp_depth_greedy_hill['depth_first_iterations'] = depth_first_iterations

df_exp_depth_greedy_hill['greedy_score'] = max(greedy_scores)
df_exp_depth_greedy_hill['greedy_iterations'] = greedy_iterations

df_exp_depth_greedy_hill['hill_climber_scores'] = max(hill_climb_scores)
df_exp_depth_greedy_hill['hill_climber_iterations'] = hill_climber_iterations

df_exp_depth_greedy_hill['sim_anneal_scores'] = max(sim_anneal_scores)
df_exp_depth_greedy_hill['sim_anneal_iterations'] = sim_anneal_iterations

df_exp_depth_greedy_hill.to_csv(path_or_buf=r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_depth_greedy_hill_try1")

# -------------------- VISUALIZE --------------------
df_exp_depth_greedy_hill = pd.read_csv(r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_depth_greedy_hill_try1")
print(df_exp_depth_greedy_hill)

iterations = [df_exp_depth_greedy_hill['depth_first_iterations'], df_exp_depth_greedy_hill['greedy_iterations'], df_exp_depth_greedy_hill['hill_climber_iterations'], df_exp_depth_greedy_hill['sim_anneal_iterations']]
scores = [df_exp_depth_greedy_hill['depth_first_score'], df_exp_depth_greedy_hill['greedy_score'], df_exp_depth_greedy_hill['hill_climber_scores'], df_exp_depth_greedy_hill['sim_anneal_scores']]

fig, ax1 = plt.subplots()

x = np.arange(5)
width = 0.40

ax1.bar(x-0.2, iterations, width)
ax1.xticks(x, ['Depth First', 'Greedy', 'Hill Climber', 'Simulated Annealing'])
ax1.xlabel("Algorithms")
ax1.ylabel("Iterations")

ax2 = ax1.twinx()
ax2.bar(x+0.2, scores, width)
ax2.ylabel("Scores")

plt.legend(["Iterations", "Scores"])
plt.title('Scores of Depth First vs Greedy vs Hill Climber vs Simulated Annealing')
plt.savefig('Depth_First_vs_Greedy_vs_Hill_Climber_vs_Simulated_Annealing_try1.pdf')
plt.show()