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
    iterations = 50

    best_score = 0
    all_scores = []
    score_dict = {}
    for split in splits:
        print(split)
        score_dict[split] = []
        for i in range(iterations):
            protein = Protein(sequence)
            greedy_test = gr.Greedy(protein, 3, splits = split)
            greedy_test.run_k()
            all_scores.append(greedy_test.protein.score)
            score_dict[split].append(greedy_test.protein.score)

            if greedy_test.protein.score < best_score:
                best_score = greedy_test.protein.score
                best_protein = greedy_test.protein


    for split in score_dict:
        sns.histplot(score_dict[split], kde=True, discrete=True, label=split)

    plt.legend(title="Splits")
    plt.title("Greedy", fontweight='bold')
    plt.xlabel("Score", loc='right')
    plt.ylabel("Frequency", loc='top')

    plt.show()

    visualize_protein(best_protein, 3)
