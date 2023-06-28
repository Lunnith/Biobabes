from ..algorithms.greedy import Greedy
from ..algorithms.simulated_annealing import SimulatedAnnealing
from ..classes.protein import Protein
from ..algorithms.hill_climber import Hill_climber
from matplotlib.lines import Line2D

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path = r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data"

# create dataframe to store results
df_exp_greedy_simanneal = pd.DataFrame()

# initialize lists to store intermediate results
split_numbers = []
score_after_greedy = []
score_after_simanneal = []

# make a dictionary with the splits as keys and the number of iterations per 'before' as the value
dict_splits = {1: 50, 2: 25, 3: 16, 4: 12, 5: 10}

# run greedy with depth with different split and before size and afterwards run simulated annealing on the folded by greedy protein
for key, value in dict_splits.items():
    for before in range(key):
        for i in range(value):
            test_protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

            greedy_protein = Greedy(test_protein, 3, splits=1, before=0)
            greedy_protein.run()

            split_numbers.append(1)
            score_after_greedy.append(greedy_protein.protein.score)

            simanneal_protein = SimulatedAnnealing(greedy_protein.protein, start_n = 10, folded = True, dimensions = 3, temperature = 10)
            simanneal_protein.run_i_iterations(greedy_protein.protein, iterations = 500, bonds = 10)
            
            score_after_simanneal.append(simanneal_protein.protein.score)

# store all results in the dataframe
df_exp_greedy_simanneal['split_numbers'] = split_numbers
df_exp_greedy_simanneal['score_after_greedy'] = score_after_greedy
df_exp_greedy_simanneal['score_after_simulated_annealing'] = score_after_simanneal

df_exp_greedy_simanneal.to_csv(path_or_buf=fr"{path}\df_exp_greedy_simanneal_bugfix")

# -------------------- VISUALIZE --------------------
df_exp_greedy_simanneal = pd.read_csv(fr"{path}\df_exp_greedy_simanneal_bugfix")
print(df_exp_greedy_simanneal)

# create a scatterplot with connected lines of the scores of the folding after greedy and the scores after simulated annealing
before_scores = df_exp_greedy_simanneal['score_after_greedy']
afters_scores = df_exp_greedy_simanneal['score_after_simulated_annealing']
split_numbs = df_exp_greedy_simanneal['split_numbers']

plt.scatter(np.zeros(len(before_scores)), before_scores, c='black')
plt.scatter(np.ones(len(afters_scores)), afters_scores, c='black')

colors = ['blue', 'red', 'green', 'orange', 'purple']

for i in range(len(before_scores)):
    plt.plot([0,1], [before_scores[i], afters_scores[i]], color=colors[split_numbs[i] - 1])

custom_lines = [Line2D([0], [0], color='blue', lw=4),
                Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='green', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='purple', lw=4)]

ax = plt.gca()
ax.invert_yaxis()
ax.legend(custom_lines, ['Split = 1', 'Split = 2', 'Split = 3', 'Split = 4', 'Split = 5'], loc='upper left')

plt.xticks([0,1], ['Score after greedy', 'Score after simulated annealing'])
plt.title('Greedy combined with Simulated Annealing')
plt.savefig(fr"{path}\Greedy_combined_with_Simulated_Annealing_BUG.pdf")
plt.show()
