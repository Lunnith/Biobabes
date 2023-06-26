from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
from ..algorithms.hill_climber import Hill_climber
from ..algorithms.simulated_annealing import SimulatedAnnealing
from ..visualization.visualize import *

sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"

# protein = Protein(sequence)
# hill_climber = Hill_climber(protein, prints=True, dimensions=3)
# folded_protein, best_score, scores, improvement = hill_climber.run_i_iterations(protein, 10, 1)
# # print(best_score)
# # visualize_protein(folded_protein, 3)
# protein = Protein(sequence)
# if folded_protein.score > -20:
#     best_protein = hill_climber.experiment(protein, 500, 2, max_n=10)
#     visualize_protein(best_protein, 3)
# else: print("Score too good:", protein.score)

protein = Protein(sequence)
sim_annealing = SimulatedAnnealing(protein, start_n=10, dimensions=3, prints=True)
folded_protein, best_score, scores, improvement = sim_annealing.run_i_iterations(protein, 100, 10)
# print(best_score)
# visualize_protein(folded_protein, 3)
# protein = Protein(sequence)
if protein.score > -20:
    best_protein, results = sim_annealing.experiment(folded_protein, 50, sample_size=2, max_n=3, sim_annealing=True, plot=False)
    print(results)
    visualize_protein(best_protein, 3)
else: print("Score too good:", protein.score)

