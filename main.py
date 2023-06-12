from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *

# initialize protein object and directions
sequence = 'HPHPHHHHHHHH'
protein = Protein(sequence)
directions_2D = set(((1, 0, 1), (-1, 0, -1), (0, 1, 2), (0, -1, -2)))

# fold protein randomly and create output
protein_folded = random_assignment(protein, directions_2D)
protein_folded.create_output()

# visualize protein
visualize_protein(protein_folded)