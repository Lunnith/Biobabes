from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.depth_first import DepthFirst
from code.algorithms.randomise import random_reassignment
#from code.algorithms.greedy import greedy
from code.algorithms import greedy as gr
import random
import copy
import matplotlib.pyplot as plt
from operator import add


if __name__ == "__main__":

    #sequence = "HCPHPCPHPCHPPCHPPCHHCHHCPPCHCPCHCPCHCHPCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHH"
    #sequence = "PCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHH"
    #sequence = "PCPHPCHCH"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCHP"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH"

    # protein = Protein(sequence)
    # greedy_test = gr.Greedy(protein, 3, splits = 3)
    # greedy_test.run_k()
    # visualize_protein(protein, 3)

    ##testing
    sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    list_scores = []
    splits = [1, 2, 3, 4, 5]
    for split in splits:
        print(split)
        for i in range(5):
            protein = Protein(sequence)
            greedy_test = gr.Greedy(protein, 3, splits = split)
            greedy_test.run_k()
            list_scores.append([split, greedy_test.protein.score])
    #visualize_protein(protein, 3)
    print(list_scores)

    dict_scores = {}
    for split in splits:
        dict_scores[split] = []
    
    for run in list_scores:
        dict_scores[run[0]].append(run[1])
    
    print(dict_scores)
