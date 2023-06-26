from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
from ..algorithms.hill_climber import Hill_climber
from ..algorithms.simulated_annealing import SimulatedAnnealing
from ..visualization.visualize import *
from ..algorithms.randomise import random_assignment
import copy

sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
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
sample_size = 10
max_n = 10
plot = False

#Now for running fast [DELETE LATER]
iterations = 10
sample_size = 2
max_n = 3

# Run the experiment for the Hill Climber
hill_climber = Hill_climber(protein=protein_hc, prints=prints, dimensions=dimensions, folded=folded)
best_protein_hc, results_hc = hill_climber.experiment(protein=protein_hc, iterations=iterations, sample_size=sample_size, max_n=max_n, plot=plot, sim_annealing=False)

# Run the experiment for the Simulated Annealing
sim_annealing = SimulatedAnnealing(protein=protein_sa, start_n=max_n, dimensions=dimensions, prints=prints, folded=folded)
best_protein_sa, results_sa = sim_annealing.experiment(protein=protein_sa, iterations=iterations, sample_size=sample_size, max_n=max_n, plot=plot, sim_annealing=True)

# Compare results
print(results_hc)
print("\n\n\n")
print(results_sa)

