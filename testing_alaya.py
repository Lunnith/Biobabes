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
import time


if __name__ == "__main__":

    #sequence = "HCPHPCPHPCHPPCHPPCHHCHHCPPCHCPCHCPCHCHPCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHHHCPCHCPCHCPPCHCCHCHCHHCPPPPCHCP"
    #sequence = "PCHCPHCPPHCPCHCPCHHPCHCPPCHCPCHCPCHPPHCHCPCHCHCHH"
    #sequence = "PCPHPCHC"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCHP"
    #sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    #sequence = "HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH"
    # sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    # protein = Protein(sequence)
    # greedy_test = gr.Greedy(protein, 3, splits = 1,  before = 0)
    # greedy_test.run()
    # visualize_protein(protein, 3)

    #testing
    #sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    #sequence = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    df_exp_greedy_2000i = pd.DataFrame()
    df_exp_greedy_2000i_total = pd.DataFrame()

    splits = []
    splits_total = []
    scores = []
    best_scores = []
    befores = []
    total_time = []

    dict_splits = {1: 2000, 2: 1000, 3: 667, 4: 500, 5: 400}
    for key, value in dict_splits.items():

        start_time = time.time()
        best_score = 1
        print(key)

        for before in range(key):
            print(before)

            for i in range(value):
                protein = Protein(sequence)
                greedy_test = gr.Greedy(protein, 3, splits = key, before = before)
                greedy_test.run()
                scores.append(greedy_test.protein.score)
                befores.append(before)
                splits.append(key)
                if greedy_test.protein.score < best_score:
                    best_protein = copy.deepcopy(greedy_test.protein)
                    best_score = best_protein.score

        end_time = time.time()

        total_time.append(end_time - start_time)
        splits_total.append(key)
        best_scores.append(best_score)

    print(best_protein)
    print(best_score)
    
    df_exp_greedy_2000i['split_numbers'] = splits
    df_exp_greedy_2000i['scores'] = scores
    df_exp_greedy_2000i['befores'] = befores

    df_exp_greedy_2000i_total['split_number'] = splits_total
    df_exp_greedy_2000i_total['best_scores'] = best_scores
    df_exp_greedy_2000i_total['total_time'] = total_time

    visualize_protein(best_protein, 3)

    print(df_exp_greedy_2000i.head())
    df_exp_greedy_2000i.to_csv(path_or_buf=r'C:\Users\alaya\Documents\Data_Greedy\df_exp_greedy_2000i')

    print(df_exp_greedy_2000i_total.head())
    df_exp_greedy_2000i_total.to_csv(path_or_buf=r'C:\Users\alaya\Documents\Data_Greedy\df_exp_greedy_2000i_total')

    # for split in score_dict:
    #     sns.histplot(score_dict[split], kde=True, discrete=True, label=split)

    # plt.legend(title="Splits")
    # plt.title("Greedy", fontweight='bold')
    # plt.xlabel("Score", loc='right')
    # plt.ylabel("Frequency", loc='top')

    # plt.show()

    # visualize_protein(best_protein, 3)
