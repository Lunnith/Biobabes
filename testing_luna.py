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
# folded_protein = random_assignment(protein, 3)
hill_climber = Hill_climber(protein, prints=False, folded=False, dimensions=2)
folded_protein = hill_climber.protein

# print("Starting score =", hill_climber.lowest_score)
folded_protein, index = hill_climber.change_bond(folded_protein)
# visualize_protein(folded_protein, 3)
folded_protein = hill_climber.refold_into_valid_state(folded_protein)

# visualize_protein(folded_protein, 3)



folded_protein, best_score, scores, improvement = hill_climber.run_i_iterations(protein, 10, 1)
# print(best_score)
# visualize_protein(folded_protein, 3)
protein = Protein(sequence)
if folded_protein.score > -20:
    best_protein = hill_climber.experiment(protein, 10, 2, max_n=10)
    visualize_protein(best_protein, 3)
else: print("Score too good:", protein.score)
print(len(best_protein.sequence_list))
print(len(sequence))