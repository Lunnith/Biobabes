from ..algorithms.greedy import Greedy
from ..algorithms.simulated_annealing import SimulatedAnnealing
from ..classes.protein import Protein
from ..algorithms.simulated_annealing import SimulatedAnnealing

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_exp_greedy_simanneal = pd.DataFrame()
split_numbers = []
score_after_greedy = []
score_after_simanneal = []
befores = []

n = 100

for i in range(n):
    for j in range(1, 6):

        for before in range(j):
            test_protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

            greedy_protein = Greedy(test_protein, 3, splits=j, before=before)
            greedy_protein.run()

            split_numbers.append(j)
            score_after_greedy.append(greedy_protein.protein.score)
            befores.append(before)

            # CHANGE TEMP BASED ON IDEAL SIMULATED ANNEALING AND ITERATIONS
            simanneal_protein = SimulatedAnnealing(greedy_protein.protein, 10, folded=True, dimensions=3, temperature=20)
            simanneal_protein.run_i_iterations(greedy_protein.protein, iterations=1000, bonds=10)
            
            score_after_simanneal.append(simanneal_protein.protein.score)

        print(f'Greedy x Simulated Annealing split {j} completed')

df_exp_greedy_simanneal['split_numbers'] = split_numbers
df_exp_greedy_simanneal['score_after_greedy'] = score_after_greedy
df_exp_greedy_simanneal['befores'] = befores
df_exp_greedy_simanneal['score_after_simulated_annealing'] = score_after_simanneal

df_exp_greedy_simanneal.to_csv(path_or_buf=r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_greedy_simanneal_try1")

# -------------------- VISUALIZE --------------------
df_exp_greedy_simanneal = pd.read_csv(r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_greedy_simanneal_try1")
print(df_exp_greedy_simanneal)

before_scores = df_exp_greedy_simanneal['score_after_greedy']
afters_scores = df_exp_greedy_simanneal['score_after_simulated_annealing']
split_numbs = df_exp_greedy_simanneal['split_numbers']

plt.scatter(np.zeros(len(before_scores)), before_scores, c='black')
plt.scatter(np.ones(len(afters_scores)), afters_scores, c='black')

colors = ['blue', 'red', 'green', 'orange', 'purple']

for i in range(len(before_scores)):
    plt.plot([0,1], [before_scores[i], afters_scores[i]], color=colors[split_numbs[i] - 1])

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='blue', lw=4),
                Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='green', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='purple', lw=4)]


ax = plt.gca()
ax.invert_yaxis()
ax.legend(custom_lines, ['Split = 1', 'Split = 2', 'Split = 3', 'Split = 4', 'Split = 5'], loc='upper left')

plt.xticks([0,1], ['Score after greedy', 'Score after simulated annealing'])
plt.title('Greedy combined with Simulated Annealing')
plt.savefig('Greedy_combined_with_Simulated_Annealing_try1.pdf')
plt.show()
