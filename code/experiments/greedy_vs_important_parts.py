from ..algorithms.greedy import Greedy
from ..algorithms.important_parts import ImportantParts
from ..classes.protein import Protein

import pandas as pd
import matplotlib.pyplot as plt

df_exp_greedy_important_parts = pd.DataFrame()

greedy_scores = []
important_parts_scores = []

# CHANGE
n = 2

for i in range(n):
    test_protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    greedy_protein = Greedy(test_protein, 3, splits=5)
    greedy_protein.run()
    greedy_scores.append(greedy_protein.protein.score)

for i in range(n):
    test_protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    important_parts_protein = ImportantParts(test_protein, 3)
    important_parts_protein.run(split_on_P=False, split_on_size=True, size=5)
    important_parts_scores.append(important_parts_protein.protein.score)

df_exp_greedy_important_parts['greedy_scores'] = greedy_scores
df_exp_greedy_important_parts['important_part_scores'] = important_parts_scores

df_exp_greedy_important_parts.to_csv(path_or_buf=r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_greedy_important_parts_try1")

# -------------------- VISUALIZE --------------------
df_exp_greedy_important_parts = pd.read_csv(r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_greedy_important_parts_try1")
print(df_exp_greedy_important_parts)

data = [df_exp_greedy_important_parts['greedy_scores'], df_exp_greedy_important_parts['important_part_scores']]

plt.boxplot(data)
plt.xlabel("Algorithm")
plt.ylabel("Score distribution")
plt.title('Greedy vs Important Parts, split size = 5, 100 iterations')

ax = plt.gca()
ax.invert_yaxis()

plt.xticks([1, 2], ['Greedy with depth', 'Important Parts'])
plt.savefig('Greedy_vs_Important Parts_try1.pdf')
plt.show()