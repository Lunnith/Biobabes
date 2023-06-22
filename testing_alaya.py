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
    sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH"

    protein = Protein(sequence)
    greedy_test = gr.Greedy(protein, 3, splits = 5,  before = 3)
    greedy_test.run()
    visualize_protein(protein, 3)

    ##testing
    sequence = "CHCCCH"
    list_scores = []
    splits = [1, 2, 3, 4, 5]
    for split in splits:
        print(split)
        for i in range(50):
            protein = Protein(sequence)
            greedy_test = gr.Greedy(protein, 3, splits = split)
            greedy_test.run_k()
            list_scores.append([split, greedy_test.protein.score])
    #visualize_protein(protein, 3)


    for split in score_dict:
        sns.histplot(score_dict[split], kde=True, discrete=True, label=split)

    plt.legend(title="Splits")
    plt.title("Greedy", fontweight='bold')
    plt.xlabel("Score", loc='right')
    plt.ylabel("Frequency", loc='top')

    plt.show()

    visualize_protein(best_protein, 3)
