from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.hill_climber import *
from code.algorithms.randomise import *
import time

#NOTE TO SELF: 2D version is buggy
# protein = Protein("HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH")
# protein = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
# sequence = "HCPHCH"



protein = Protein(sequence)
hill_climber = Hill_climber(protein, 3)
folded_protein = hill_climber.protein

# print("Starting score =", hill_climber.lowest_score)
folded_protein, index = hill_climber.change_bond(protein)
folded_protein = hill_climber.refold(protein, index)
# visualize_protein(folded_protein, 3)
# folded_protein = hill_climber.refold_into_valid_state(protein)

# visualize_protein(folded_protein, 3)



folded_protein, best_score, scores, improvement = hill_climber.run_n_iterations(protein, 10, 1)
# print(best_score)
# visualize_protein(folded_protein, 3)
if folded_protein.score > -20:
    best_protein = hill_climber.experiment(folded_protein, 5, 3, max_n=5)
    visualize_protein(best_protein, 3)
else: print("Score too good:", folded_protein.score)

