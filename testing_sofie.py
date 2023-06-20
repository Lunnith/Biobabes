from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *
from code.algorithms.depth_first import DepthFirst
from code.algorithms.important_parts import ImportantParts

import matplotlib.pyplot as plt
import time

if __name__ == "__main__":

    # initialize protein object and directions
    sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    test_protein = Protein(sequence)

    #fold protein randomly and create output
    # depth_first = DepthFirst(test_protein, 3)
    # depth_first.run(P_pruning=True, directions_pruning=True)
    # depth_first.protein.create_output()
    # visualize_protein(depth_first.protein, 3)

    # number_of_checked_states = depth_first.number_of_states_checked
    # list_of_best_scores = depth_first.best_scores_list

    # plt.plot(number_of_checked_states, list_of_best_scores)
    # plt.show()

    start = time.time()

    important_parts = ImportantParts(test_protein, 3)
    important_parts.run_in_parts(n=1000, split_on_P=False, split_on_size=True, size=10)
    important_parts.protein.create_output()

    end = time.time()
    print(f'Time: {end - start}')
    
    visualize_protein(important_parts.protein, 3)