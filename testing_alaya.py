from code.classes.protein import Protein
from code.visualization.visualize import *
#from code.algorithms.greedy import greedy
from code.algorithms import greedy as gr
import random
import copy
import matplotlib.pyplot as plt
from operator import add


if __name__ == "__main__":

    sequence = "HCPHPCPHPCHPPCHPPCHHCHHCPPCHCPCHCPCHCHPCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHH"
    #sequence = "HPCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHH"
    #sequence = "HCPHPCPHPCHCH"
    protein = Protein(sequence)
    greedy_test = gr.Greedy(protein, 3)
    greedy_test.run()
    visualize_protein(protein, 3)

