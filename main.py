from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *

# Score en bonds gaan nog niet goed bij uitvoeren van de random_assignment functie
# mogelijke oplossing is random assignment functie een protein object aan laten maken (werkt bij make_and_fold_protein functie)

# initialize protein object and directions
sequence = 'HPHPHHHHHHHH'
protein = Protein(sequence)
directions_2D = set(((1, 0, 1), (-1, 0, -1), (0, 1, 2), (0, -1, -2)))

# fold protein randomly and create output
protein = random_assignment(protein, directions_2D)
protein.create_output()

# visualize protein
visualize_protein(protein)