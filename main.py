from code.classes.protein import Protein

from code.algorithms.randomise import *
from code.algorithms.greedy import Greedy
from code.algorithms.depth_first import DepthFirst
from code.algorithms.important_parts import ImportantParts
from code.algorithms.hill_climber import Hill_climber
from code.algorithms.simulated_annealing import SimulatedAnnealing

from code.visualization.visualize import *

if __name__ == "__main__":

    # all proteins given in the case
    test_protein_A = Protein("HHPHHHPH")
    test_protein_B = Protein("HHPHHHPHPHHHPH")
    test_protein_C = Protein("HPHPPHHPHPPHPHHPPHPH")
    test_protein_D = Protein("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    test_protein_E = Protein("HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH")
    test_protein_F = Protein("PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP")
    test_protein_G = Protein("CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC")
    test_protein_H = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
    test_protein_I = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

    # ---------------------------- Random ----------------------------
    # random_protein, score_list = random_reassignment(test_protein_A, 3, k=100000)
    
    # print(f'Value of the folding after Randomized Assignment:'
    #       f'{random_protein.score}')

    # ---------------------------- Greedy ----------------------------
#     greedy_protein = Greedy(test_protein_A, 3, splits=1)
#     greedy_protein.run()

#     print(f'Value of the folding after Greedy:'
#           f'{greedy_protein.protein.score}')
    
    # ---------------------------- Greedy with beam search ----------------------------
#     greedy_beam_protein = Greedy(test_protein_A, 3, splits=3)
#     greedy_beam_protein.run()

#     print(f'Value of the folding after Greedy with beam search:'
#           f'{greedy_beam_protein.protein.score}')

    # ---------------------------- Depth First with P and directions pruning ----------------------------
    # Note: for longer proteins (+16 aminoacids) this will take longer than an hour until forever.

#     depth_first_protein = DepthFirst(test_protein_A, 3)
#     depth_first_protein.run(P_pruning=True, directions_pruning=True)

#     print(f'Value of the folding after Greedy:'
#           f'{depth_first_protein.protein.score}')
    
    # ---------------------------- Depth First - Important Parts ----------------------------
#     important_parts_protein = ImportantParts(test_protein_A, 3)
#     important_parts_protein.run(iterations=1000, split_on_P=True, split_on_size=False)

#     print(f'Value of the folding after Greedy:'
#           f'{important_parts_protein.protein.score}')
    
    # ---------------------------- Hill Climber ----------------------------
#     hill_climber_protein = Hill_climber(test_protein_A, 3)
#     hill_climber_protein.run_i_iterations(test_protein_A, iterations=1000, bonds=1)

#     print(f'Value of the folding after Hill Climber:'
#           f'{hill_climber_protein.protein.score}')
    
    # ---------------------------- Simulated Annealing ----------------------------
    simanneal_protein = SimulatedAnnealing(test_protein_I, 10, 3,  temperature=20)
    simanneal_protein.run_i_iterations(test_protein_I, iterations=1000, bonds=10)

    print(f'Value of the folding after Simulated Annealing:'
          f'{simanneal_protein.protein.score}')
    
    # ---------------------------- Visualization ----------------------------
    visualize_protein(simanneal_protein.protein, 3)

