from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.hill_climber import *

protein = Protein("HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH")
# protein = Protein("HPPCCHPPH")
hill_climber = Hill_climber(protein)
folded_protein = hill_climber.protein
print("Starting score =", hill_climber.lowest_score)
# hill_climber.change_one_bond(protein)
visualize_protein(folded_protein, 3)



folded_protein, best_score = hill_climber.run_n_iterations(protein, 500, 1)
print(best_score)
visualize_protein(folded_protein, 3)

