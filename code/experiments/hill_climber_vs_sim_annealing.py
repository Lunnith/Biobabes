from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
from ..algorithms.hill_climber import Hill_climber
from ..algorithms.simulated_annealing import SimulatedAnnealing
from ..visualization.visualize import *
from ..algorithms.randomise import random_assignment
import copy

sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
# note: You can now run hill_climber.experiment() and simulated_annealing.experiment()
#       And it will return a dataframe with the lines that would be plotted
#   Now, create a function to compare these.

#Prepare protein
protein = random_assignment(Protein(sequence), 3)
while len(sequence) != len(protein.sequence_list) or protein.score < -10:
    protein = random_assignment(Protein(sequence), 3)
protein_hc = copy.deepcopy(protein)
protein_sa = copy.deepcopy(protein)

print(f"Starting score = {protein.score}")

# Make sure to use the same kwargs for both experiments
dimensions = 3
prints = True
folded = True

iterations = 500
sample_size = 50
max_n = 10
plot = False

# Run the experiment for the Hill Climber
hill_climber = Hill_climber(protein=protein_hc, prints=prints, dimensions=dimensions, folded=folded)
best_protein_hc, results_hc = hill_climber.experiment(protein=protein_hc, iterations=iterations, sample_size=sample_size, max_n=max_n, plot=plot, sim_annealing=False)

# Run the experiment for the Simulated Annealing
sim_annealing = SimulatedAnnealing(protein=protein_sa, start_n=max_n, dimensions=dimensions, prints=prints, folded=folded)
best_protein_sa, results_sa = sim_annealing.experiment(protein=protein_sa, iterations=iterations, sample_size=sample_size, max_n=max_n, plot=plot, sim_annealing=True)

# unpack results
# hc_index = []
# for n in range(1, max_n+1):
#     hc_index.append(f"n={n}")
# results_hc.index = hc_index
results_hc.loc["Algorithm"] = 'Hill Climber'
df_hc = copy.deepcopy(results_hc.T[0:-1])
best_hc = copy.deepcopy(results_hc.T.iloc[-1])

# sa_index = []
# for n in range(0, max_n):
#     sa_index.append(f"n={n}")
# results_sa.index = sa_index
results_sa.loc["Algorithm"] = 'Simulated Annealing'
df_sa = copy.deepcopy(results_sa.T[0:-1])
best_sa = copy.deepcopy(results_sa.T.iloc[-1])

df_both_messy = pd.concat((df_hc, df_sa))
# print(df_both_messy)
df_both = df_both_messy.melt(id_vars=['Algorithm'], var_name='Bonds changed', value_name='Score')
# df_both.index = df_both['n_used']
# df_both = df_both.drop(labels='n_used', axis=1)
print(df_both)

#Plot results
# fig, axs = plt.subplots(1, 2, figsize=(12, 5))
# ax = 0
used_n = range(0, max_n)

ax = sns.boxplot(data=df_both, y='Score', x='Bonds changed', hue='Algorithm')
ax.invert_yaxis()
plt.show()
