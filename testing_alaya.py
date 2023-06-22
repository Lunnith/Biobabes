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

    sequence = "HCPHPCPHPCHPPCHPPCHHCHHCPPCHCPCHCPCHCHPCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHHHCPCHCPCHCPPCHCCHCHCHHCPPPPCHCP"
    #sequence = "PCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHH"
    #sequence = "PCPHPCHC"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCHP"
    #sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH"

    # protein = Protein(sequence)
    # greedy_test = gr.Greedy(protein, 2, splits = 2,  before = 0)
    # greedy_test.run()
    # visualize_protein(protein, 2)

    ##testing
    #sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    #sequence = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    best_score = 1
    list_scores = []
    splits = [1, 2, 3, 4, 5]
    for split in splits:
        for before in range(split):
            if split == 1 or split == 2 or split == 3:
                extra_iterations = 100
            else:
                extra_iterations = 0
                for i in range(50//split + extra_iterations):
                    protein = Protein(sequence)
                    greedy_test = gr.Greedy(protein, 3, splits = split, before = before)
                    greedy_test.run()
                    list_scores.append([split, greedy_test.amino_before, greedy_test.protein.score])
                    if greedy_test.protein.score < best_score:
                        best_protein = copy.deepcopy(greedy_test.protein)
                        best_score = best_protein.score
    print(best_protein)
    print(best_score)
    print(list_scores)

    visualize_protein(best_protein, 3)


    # for split in score_dict:
    #     sns.histplot(score_dict[split], kde=True, discrete=True, label=split)

    # plt.legend(title="Splits")
    # plt.title("Greedy", fontweight='bold')
    # plt.xlabel("Score", loc='right')
    # plt.ylabel("Frequency", loc='top')

    # plt.show()

    # visualize_protein(best_protein, 3)
