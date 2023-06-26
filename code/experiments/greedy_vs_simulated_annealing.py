from ..algorithms.greedy import Greedy
from ..algorithms.simulated_annealing import SimulatedAnnealing

import pandas as pd
import matplotlib.pyplot as plt

df_exp_greedy_simanneal = pd.DataFrame()
split_numbers = []
score_after_greedy = []
score_after_simanneal = []
befores = []

# CHANGE
n = 1

for i in range(1, 5):
    for j in range(n):

        before_list = []

        for number in range(i):
            if number < i:
                before_list.append(number)

        for before in before_list:
            test_protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

            greedy_protein = Greedy(test_protein, 3, splits=i, before=before)
            greedy_protein.run()

            print(f'Value of the folding after Greedy:'
                f'{greedy_protein.protein.score}')

            split_numbers.append(i)
            score_after_greedy.append(greedy_protein.protein.score)
            befores.append(before)

            # CHANGE TEMP BASED ON IDEAL SIMULATED ANNEALING
            simanneal_protein = SimulatedAnnealing(greedy_protein.protein, 10, folded=True, dimensions=3, temperature=20)
            simanneal_protein.run_i_iterations(greedy_protein.protein, iterations=10000, bonds=10)

            print(f'Value of the folding after Simulated Annealing x Greedy:'
                    f'{simanneal_protein.protein.score}')
            
            score_after_simanneal.append(simanneal_protein.protein.score)

df_exp_greedy_simanneal['split_numbers'] = split_numbers
df_exp_greedy_simanneal['score_after_greedy'] = score_after_greedy
df_exp_greedy_simanneal['befores'] = befores
df_exp_greedy_simanneal['score_after_simulated_annealing'] = score_after_simanneal

print(df_exp_greedy_simanneal.head())
df_exp_greedy_simanneal.to_csv(path_or_buf=r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_greedy_simanneal_try1")

# -------------------- VISUALIZE --------------------
df_exp_greedy_simanneal = pd.read_csv(r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_greedy_simanneal_try1")
print(df_exp_greedy_simanneal)

befores = df_exp_greedy_simanneal['scores_after_greedy']
afters = df_exp_greedy_simanneal['scores_after_simanneal']

plt.scatter(np.zeros(len(befores)), befores)
plt.scatter(np.ones(len(afters)), afters)

for i in range(len(befores)):
    plt.plot( [0,1], [befores[i], afters[i]])
    plt.legend(f'Split: {df_exp_greedy_simanneal['split_numbers'][i]} and before: {df_exp_greedy_simanneal['befores'][i]}')

plt.xticks([0,1], ['Score after greedy', 'Score after simulated annealing'])
plt.title('Greedy combined with Simulated Annealing')
plt.savefig('Greedy_combined_with_Simulated_Annealing_try1.pdf')
plt.show()
