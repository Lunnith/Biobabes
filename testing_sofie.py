from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *
from code.algorithms.depth_first import DepthFirst
from code.algorithms.depth_first import RecognizePatterns

import matplotlib.pyplot as plt

if __name__ == "__main__":

    # initialize protein object and directions
    sequence = "CCHHPPPCHPCPCCCHHHPCH"
    test_protein = Protein(sequence)

    # fold protein randomly and create output
    # depth_first = DepthFirst(test_protein, 3)
    # depth_first.run(P_pruning=False, directions_pruning=True)
    # depth_first.protein.create_output()
    # visualize_protein(depth_first.protein, 3)

    # number_of_checked_states = depth_first.number_of_states_checked
    # list_of_best_scores = depth_first.best_scores_list

    # plt.plot(number_of_checked_states, list_of_best_scores)
    # plt.show()


    recognize_patterns = RecognizePatterns(test_protein, 3)
    recognize_patterns.run_in_parts()
    recognize_patterns.protein.create_output()
    visualize_protein(recognize_patterns.protein, 3)