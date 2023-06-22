from ..algorithms.greedy import Greedy
from ..algorithms.hill_climber import Hill_climber

from ..visualization import *

if __name__ == "__main__":   

    test_protein_I = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    
    greedy_protein = Greedy(test_protein_I, 3, splits=1)
    greedy_protein.run()

    print(f'Value of the folding after Greedy:'
          f'{greedy_protein.protein.score}')
    

    hill_climber_protein = Hill_climber(greedy_protein, 3)
    hill_climber_protein.run_i_iterations(greedy_protein, iterations=1000, bonds=1)

    print(f'Value of the folding after Hill Climber x Greedy:'
          f'{hill_climber_protein.protein.score}')
    
    visualize_protein(hill_climber_protein.protein, 3)