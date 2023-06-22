from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.depth_first import DepthFirst
from code.algorithms.randomise import random_reassignment
#from code.algorithms.greedy import greedy
from code.algorithms import greedy as gr
import random
import copy
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from operator import add


if __name__ == "__main__":

    #sequence = "HCPHPCPHPCHPPCHPPCHHCHHCPPCHCPCHCPCHCHPCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHHHCPCHCPCHCPPCHCCHCHCHHCPPPPCHCP"
    #sequence = "PCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHH"
    #sequence = "PCPHPCHC"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCHP"
    #sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH"
    #sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    # protein = Protein(sequence)
    # greedy_test = gr.Greedy(protein, 3, splits = 4,  before = 0)
    # greedy_test.run()
    # visualize_protein(protein, 3)

    #testing
    #sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    #sequence = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    best_score = 1
    list_best_scores = []
    list_scores = []
    dict_splits = {1: 10000}
                   #, 2: 1000, 3: 300, 4: 20, 5: 10}
    for key, value in dict_splits.items():
        for before in range(key):
            for i in range(value):
                protein = Protein(sequence)
                greedy_test = gr.Greedy(protein, 3, splits = key, before = before)
                greedy_test.run()
                list_scores.append([key, greedy_test.protein.score])
                if greedy_test.protein.score < best_score:
                    best_protein = copy.deepcopy(greedy_test.protein)
                    best_score = best_protein.score
        list_best_scores.append([key, best_protein.score])
    print(best_protein)
    print(best_score)
    print(list_scores)
    print(list_best_scores)

    visualize_protein(best_protein, 3)


    # for split in score_dict:
    #     sns.histplot(score_dict[split], kde=True, discrete=True, label=split)

    # plt.legend(title="Splits")
    # plt.title("Greedy", fontweight='bold')
    # plt.xlabel("Score", loc='right')
    # plt.ylabel("Frequency", loc='top')

    # plt.show()

    # visualize_protein(best_protein, 3)
