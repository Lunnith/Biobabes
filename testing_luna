from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.hill_climber import *

protein = Protein("HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH")
hill_climber = Hill_climber(protein)
folded_protein = hill_climber.protein
visualize_protein(folded_protein, 3)



folded_protein, best_score = hill_climber.run_n_iterations(protein, 5, 3)
print(best_score)
visualize_protein(folded_protein, 3)


