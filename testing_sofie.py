from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *
from code.algorithms.depth_first import DepthFirst
from code.algorithms.important_parts import ImportantParts
from code.algorithms.simulated_annealing import SimulatedAnnealing

import matplotlib.pyplot as plt
import time

if __name__ == "__main__":

    # initialize protein object and directions
    sequence = "HCPHP"
    test_protein = Protein(sequence)
    depth_first = DepthFirst(test_protein, 3)
    depth_first.run(P_pruning=False, directions_pruning=False)
    depth_first.protein.create_output()
    visualize_protein(depth_first.protein, 3)

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

    # important_parts = ImportantParts(test_protein, 3)
    # important_parts.run()
    # important_parts.protein.create_output()

    # # end = time.time()
    # # print(f'Time: {end - start}')
    
    # visualize_protein(important_parts.protein, 3)

    # simanneal_protein = SimulatedAnnealing(test_protein, 3, temperature=5)
    # simanneal_protein.run_n_iterations(test_protein, iterations=5000, bonds=1)

    # visualize_protein(simanneal_protein.protein, 3)
    # print(f'Value of the folding after Simulated Annealing:'
    #       f'{simanneal_protein.protein.score}')