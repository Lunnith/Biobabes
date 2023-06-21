from code.classes.protein import Protein

from code.algorithms.randomise import *
from code.algorithms.greedy import Greedy
from code.algorithms.depth_first import DepthFirst
from code.algorithms.important_parts import ImportantParts
from code.algorithms.hill_climber import Hill_climber
from code.algorithms.simulated_annealing import SimulatedAnnealing

from code.visualization.visualize import *

if __name__ == "__main__":

    # create a protein with a given sequence
    sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    test_protein = Protein(sequence)

    # ---------------------------- Random ----------------------------
    # random_protein, score_list = random_reassignment(test_protein, 3, k=100000)
    
    # print(f'Value of the folding after Randomized Assignment:'
    #       f'{random_protein.protein.score}')

    # ---------------------------- Greedy ----------------------------
    # greedy_protein = Greedy(test_protein, 3, splits=1)
    # greedy_protein.run()

    # print(f'Value of the folding after Greedy:'
    #       f'{greedy_protein.protein.score}')
    
    # ---------------------------- Greedy with beam search ----------------------------
    # greedy_beam_protein = Greedy(test_protein, 3, splits=9)
    # greedy_beam_protein.run()

    # print(f'Value of the folding after Greedy with beam search:'
    #       f'{greedy_beam_protein.protein.score}')

    # ---------------------------- Depth First with P and directions pruning ----------------------------
    # Note: for longer proteins (+16 aminoacids) this will take longer than an hour until forever.
    # depth_first_protein = DepthFirst(test_protein, 3)
    # depth_first_protein.run(P_pruning=True, directions_pruning=True)

    # print(f'Value of the folding after Greedy:'
    #       f'{depth_first_protein.protein.score}')
    
    # ---------------------------- Depth First - Important Parts ----------------------------
    # important_parts_protein = ImportantParts(test_protein, 3)
    # important_parts_protein.run_in_parts(n=1000, split_on_P=True, split_on_size=False)

    # print(f'Value of the folding after Greedy:'
    #       f'{important_parts_protein.protein.score}')
    
    # ---------------------------- Hill Climber ----------------------------
    # hill_climber_protein = Hill_climber(test_protein, 3)
    # hill_climber_protein.run_n_iterations(test_protein, iterations=1000, bonds=1)

    # print(f'Value of the folding after Hill Climber:'
    #       f'{hill_climber_protein.protein.score}')
    
    # ---------------------------- Simulate Annealing ----------------------------
    simanneal_protein = SimulatedAnnealing(test_protein, 3, temperature=5)
    simanneal_protein.run_n_iterations(test_protein, iterations=1000, bonds=1)

    print(f'Value of the folding after Simulated Annealing:'
          f'{simanneal_protein.protein.score}')
    
    # ---------------------------- Visualization ----------------------------
    visualize_protein(simanneal_protein.protein, 3)

