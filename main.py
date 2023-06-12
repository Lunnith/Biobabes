from code.classes.protein import Protein
from code.classes.aminoacid import Aminoacid
from code.visualization.visualize import *
from code.algorithms.randomise import *

if __name__ == "__main__":

    # initialize protein object and directions
    sequence = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'
    protein = Protein(sequence)

    # fold protein randomly and create output
    folded_protein, score_list = random_reassignment(protein, 3, 100)
    folded_protein.create_output()

    # visualize protein
    visualize_scores(score_list)
    visualize_protein(folded_protein, 3)