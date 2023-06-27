from ..algorithms.greedy import Greedy
from ..algorithms.important_parts import ImportantParts
from ..classes.protein import Protein

import pandas as pd
import matplotlib.pyplot as plt

path = r"C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data"

# create dataframe to store results
df_exp_greedy_important_parts = pd.DataFrame()

# initialize lists
greedy_scores = []
important_parts_scores = []

# run important parts algorithm with split size 5 n iterations
n = 1000

for i in range(n):
    test_protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    important_parts_protein = ImportantParts(test_protein, 3)
    important_parts_protein.run(iterations=1000, split_on_P=False, split_on_size=True, size=5)
    important_parts_scores.append(important_parts_protein.protein.score)

print(f'Important Parts done!')

df_exp_greedy_important_parts.to_csv(path_or_buf=fr"{path}\df_exp_greedy_important_parts_complete")

# run greedy with depth algorithm with split size 5 n iterations
n = 500

for i in range(n):
    for before in range(5):
        test_protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
        greedy_protein = Greedy(test_protein, 3, splits=5, before=before)
        greedy_protein.run()
        greedy_scores.append(greedy_protein.protein.score)

print('Greedy done!')

# store lists in dataframe and save dataframe
df_exp_greedy_important_parts['greedy_scores'] = greedy_scores
df_exp_greedy_important_parts['important_part_scores'] = important_parts_scores

df_exp_greedy_important_parts.to_csv(path_or_buf=fr"{path}\df_exp_greedy_important_parts_complete")

# -------------------- VISUALIZE --------------------
df_exp_greedy_important_parts = pd.read_csv(fr"{path}\df_exp_greedy_important_parts_complete")

# plot score distributions of both algorithms in a boxplot
data = [df_exp_greedy_important_parts['greedy_scores'], df_exp_greedy_important_parts['important_part_scores']]

plt.boxplot(data)
plt.xlabel("Algorithm")
plt.ylabel("Score distribution")
plt.title('Greedy vs Important Parts, split size = 5, 1000 iterations')

ax = plt.gca()
ax.invert_yaxis()

plt.xticks([1, 2], ['Greedy with depth', 'Important Parts'])
plt.savefig(fr"{path}\Greedy_vs_Important_Parts_complete.pdf")
plt.show()