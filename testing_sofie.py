from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.randomise import *
from code.algorithms.depth_first import DepthFirst
from code.algorithms.breadth_first import BreadthFirst

if __name__ == "__main__":

    # initialize protein object and directions
    sequence = "CHPHCHHC"
    test_protein = Protein(sequence)

    # fold protein randomly and create output
    depth_first = DepthFirst(test_protein, 3)
    depth_first.run(P_pruning=True, directions_pruning=True)
    depth_first.protein.create_output()
    visualize_protein(depth_first.protein, 3)