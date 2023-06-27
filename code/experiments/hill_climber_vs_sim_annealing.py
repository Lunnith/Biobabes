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
best_protein_hc, results_hc, best_results_hc = hill_climber.experiment(protein=protein_hc, iterations=iterations, sample_size=sample_size, max_n=max_n, plot=plot, sim_annealing=False, result_each_sample=True)

# Run the experiment for the Simulated Annealing
sim_annealing = SimulatedAnnealing(protein=protein_sa, start_n=max_n, dimensions=dimensions, prints=prints, folded=folded)
best_protein_sa, results_sa, best_results_sa = sim_annealing.experiment(protein=protein_sa, iterations=iterations, sample_size=sample_size, max_n=max_n, plot=plot, sim_annealing=True, result_each_sample=True)

# unpack results per sample into one big dataframe of results
unpacked_results_hc = {}
for n in best_results_hc:
    unpacked_results_hc[n] = []
    for sample in best_results_hc[n]:
        unpacked_results_hc[n].append(best_results_hc[n][sample])
df_results_hc = pd.DataFrame.from_dict(unpacked_results_hc).T

unpacked_results_sa = {}
for n in best_results_sa:
    unpacked_results_sa[n] = []
    for sample in best_results_sa[n]:
        unpacked_results_sa[n].append(best_results_sa[n][sample])
df_results_sa = pd.DataFrame.from_dict(unpacked_results_sa).T

#Combine results
df_results_hc.loc["Algorithm"] = 'Hill Climber'
df_hc = copy.deepcopy(df_results_hc.T[0:-1])
best_hc = copy.deepcopy(df_results_hc.T.iloc[-1])

df_results_sa.loc["Algorithm"] = 'Simulated Annealing'
df_sa = copy.deepcopy(df_results_sa.T[0:-1])
best_sa = copy.deepcopy(df_results_sa.T.iloc[-1])

df_both_messy = pd.concat((df_hc, df_sa))
df_both_messy = pd.DataFrame(df_both_messy)
df_both = df_both_messy.melt(id_vars=['Algorithm'], var_name='Bonds changed', value_name='Score')
print(df_both)

#Plot results
used_n = range(0, max_n)

ax = sns.boxplot(data=df_both, y='Score', x='Bonds changed', hue='Algorithm')
ax.invert_yaxis()
plt.show()
