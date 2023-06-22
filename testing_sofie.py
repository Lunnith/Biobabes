from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *
from code.algorithms.depth_first import DepthFirst
from code.algorithms.important_parts import ImportantParts
from code.algorithms.simulated_annealing import SimulatedAnnealing
from code.algorithms.greedy import Greedy
from code.algorithms.hill_climber import Hill_climber

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

# ------------------------------------- EXPERIMENT 1 ----------------------------------------------
# EXPERIMENT depth first score distribution with and without pruning
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


# ------------------------------------- EXPERIMENT 2 ----------------------------------------------
# EXPERIMENT greedy x simulated annealing

# df_exp_greedy_simanneal = pd.DataFrame()
# split_numbers = []
# score_after_greedy = []
# score_after_simanneal = []
# befores = []

# for i in range(1, 5):
#     for j in range(10):
#         test_protein_I = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

#         before_list = []

#         for number in range(i):
#             if number < i:
#                 before_list.append(number)
        
#         print(before_list)

#         for before in before_list:
#             greedy_protein = Greedy(test_protein_I, 3, splits=i, before)
#             greedy_protein.run()

#             print(f'Value of the folding after Greedy:'
#                 f'{greedy_protein.protein.score}')

#             split_numbers.append(i)
#             score_after_greedy.append(greedy_protein.protein.score)
#             befores.append(before)

#             simanneal_protein = SimulatedAnnealing(greedy_protein.protein, 10, folded=True, dimensions=3, temperature=20)
#             simanneal_protein.run_i_iterations(greedy_protein.protein, iterations=10000, bonds=10)

#             print(f'Value of the folding after Simulated Annealing x Greedy:'
#                     f'{simanneal_protein.protein.score}')
            
#             score_after_simanneal.append(simanneal_protein.protein.score)

# df_exp_greedy_simanneal['split_numbers'] = split_numbers
# df_exp_greedy_simanneal['score_after_greedy'] = score_after_greedy
# df_exp_greedy_simanneal['befores'] = befores
# df_exp_greedy_simanneal['score_after_simulated_annealing'] = score_after_simanneal

# print(df_exp_greedy_simanneal.head())
# df_exp_greedy_simanneal.to_csv(path_or_buf=r'C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_greedy_simanneal_complete')



# ------------------------------------- EXPERIMENT 3 ----------------------------------------------
# EXPERIMENT: depth first vs greedy/hill climb in 2d
# df_exp_depth_greedy_hill = pd.DataFrame()

# greedy_scores = []
# hill_climb_scores = []
# sim_anneal_scores = []

# greedy_iterations = 0
# hill_climber_iterations = 0
# depth_first_iterations = 1
# sim_anneal_iterations = 0

# sequence = "HHPHHHPHPHHHPH"
# test_protein = Protein(sequence)

# depth_first = DepthFirst(test_protein, 2)
# depth_first.run(P_pruning=True, directions_pruning=True)

# # CHANGE
# n = 0

# for i in range(n):
#     for i in range(1, 7):

#         before_list = []

#         for number in range(i):
#             if number < i:
#                 before_list.append(number)

#         for before in before_list:
#             greedy_protein = Greedy(test_protein_I, 2, splits=i, before)
#             greedy_protein.run()

#             greedy_scores.append(greedy_protein.protein.score)
#             greedy_iterations += 1
    
# for i in range(n):
#     hill_climber_protein = Hill_climber(test_protein_A, 2)
#     hill_climber_protein.run_i_iterations(test_protein_A, iterations=500, bonds=1)

#     hill_climber_scores.append(hill_climber_protein.protein.score)
#     hill_climber_iterations += 1

# for i in range(n):
#     for i in range(1, 10):
#         simanneal_protein = SimulatedAnnealing(test_protein_I, i, 2, temperature=20)
#         simanneal_protein.run_i_iterations(test_protein_I, iterations=500, bonds=i)

#         sim_anneal_scores.append(simanneal_protein.protein.score)
#         sim_anneal_iterations += 1


# df_exp_depth_greedy_hill = pd.DataFrame()

# df_exp_depth_greedy_hill['depth_first_score'] = depth_first.protein.score
# df_exp_depth_greedy_hill['depth_first_iterations'] = depth_first_iterations

# df_exp_depth_greedy_hill['greedy_score'] = max(greedy_scores)
# df_exp_depth_greedy_hill['greedy_iterations'] = greedy_iterations

# df_exp_depth_greedy_hill['hill_climber_scores'] = max(hill_climb_scores)
# df_exp_depth_greedy_hill['hill_climber_iterations'] = hill_climber_iterations

# df_exp_depth_greedy_hill['sim_anneal_scores'] = max(sim_anneal_scores)
# df_exp_depth_greedy_hill['sim_anneal_iterations'] = sim_anneal_iterations

# df_exp_greedy_simanneal.to_csv(path_or_buf=r'C:\Users\sofie\minorAI\Algoritmen en Heuristieken\data\df_exp_greedy_simanneal_complete')


# ------------------------------------- EXPERIMENT 4 ----------------------------------------------
# EXPERIMENT greedy vs important parts


    # greedy_protein = Greedy(test_protein_I, 3, splits=3)
    # greedy_protein.run()

    # important_parts = ImportantParts(test_protein, 3)
    # important_parts.run()
    # important_parts.protein.create_output()


# ------------------------------------- VISUALIZE ----------------------------------------------