from ..algorithms.randomise import random_reassignment
import matplotlib.pyplot as plt
from ..visualization.visualize import *
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid


sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
protein = Protein(sequence)
best_protein, scores, best_score = random_reassignment(protein, 3, k = 100000)

print(scores)
print('best score', best_score)
visualize_protein(best_protein, 3)
