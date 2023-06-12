from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *

if __name__ == "__main__":

    # initialize protein object and directions
    sequence = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'
    protein = Protein(sequence)
    directions_2D = set(((1, 0, 1), (-1, 0, -1), (0, 1, 2), (0, -1, -2)))

    # fold protein randomly and create output
    folded_protein, score_list = random_reassignment(protein, directions_2D, 100)
    folded_protein.create_output()

    # visualize protein
    visualize_scores(score_list)
    visualize_protein(folded_protein)