from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *
from code.algorithms.depth_first import DepthFirst
from code.algorithms.important_parts import ImportantParts

import matplotlib.pyplot as plt
import time

if __name__ == "__main__":

    # initialize protein object and directions
    sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    test_protein = Protein(sequence)
    # depth_first = DepthFirst(test_protein, 3)
    # depth_first.run(P_pruning=False, directions_pruning=False)
    # visualize_protein(depth_first.protein, 3)

    # experiment with and without pruning
    # scores_wo_pruning = []
    # scores_p_pruning = []
    # scores_direction_pruning = []

    # for i in range(5):
    #     test_protein = Protein(sequence)
    #     depth_first = DepthFirst(test_protein, 3)
    #     depth_first.run(P_pruning=False, directions_pruning=False)
    #     scores_wo_pruning.append(depth_first.protein.score)

    #     print('Depth First run completed')

    # print('Without pruning completed')

    # for i in range(5):
    #     test_protein = Protein(sequence)
    #     depth_first = DepthFirst(test_protein, 3)
    #     depth_first.run(P_pruning=True, directions_pruning=False)
    #     scores_p_pruning.append(depth_first.protein.score)
    #     print('Depth First run completed')

    # print('With p-pruning completed')

    # for i in range(5):
    #     test_protein = Protein(sequence)
    #     depth_first = DepthFirst(test_protein, 3)
    #     depth_first.run(P_pruning=False, directions_pruning=True)
    #     scores_direction_pruning.append(depth_first.protein.score)
    #     print('Depth First run completed')

    # print('With directions pruning completed')

    # color_list = ['r', 'b', 'g']

    # plt.hist(x=[scores_wo_pruning, scores_p_pruning, scores_direction_pruning], color=color_list, bins=3)
    # plt.xlabel('Score')
    # plt.ylabel('Frequency')
    # plt.legend(["Without pruning", "P-pruning", "Direction pruning"])    
    # plt.title('Depth First with different pruning types', fontweight="bold")
    # plt.xlim(right=0)
    # plt.show()

    # visualize_protein(depth_first.protein, 3)

    # number_of_checked_states = depth_first.number_of_states_checked
    # list_of_best_scores = depth_first.best_scores_list

    # plt.plot(number_of_checked_states, list_of_best_scores)
    # plt.show()

    # start = time.time()

    important_parts = ImportantParts(test_protein, 3)
    important_parts.run_in_parts(n=1000, split_on_P=True, split_on_size=False)
    important_parts.protein.create_output()

    # end = time.time()
    # print(f'Time: {end - start}')
    
    visualize_protein(important_parts.protein, 3)