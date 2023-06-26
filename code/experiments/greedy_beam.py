from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms import greedy as gr
import random
import copy
import pandas as pd
from operator import add
import time    


sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
df_exp_greedy_big_parts = pd.DataFrame()
df_exp_greedy_big_parts_total = pd.DataFrame()

# make empty lists to fill during the experiment
splits = []
splits_total = []
scores = []
best_scores = []
befores = []
total_time = []

# make a dictionary with the splits as keys and the number of iterations per 'before' as the value
dict_splits = {1: 2000, 2: 1000, 3: 667, 4: 500, 5: 400}

for key, value in dict_splits.items():
    
    # start the clock for measuring the time of one every split
    start_time = time.time()

    # set best score to 1
    best_score = 1

    # print the key to see where you are in the experiment
    print(key)

    # for every before in the range of the split do the number of iterations for that split
    for before in range(key):
        print(before)

        for i in range(value):
            print(i)

            # make a protein, initialize greedy with the right parameters and run the algorithm
            protein = Protein(sequence)
            greedy_test = gr.Greedy(protein, 3, splits = key, before = before)
            greedy_test.run()

            # append for each iterations the score, the before and the split
            scores.append(greedy_test.protein.score)
            befores.append(before)
            splits.append(key)

            # if the score is lower than the best score, set best score at protein score
            if greedy_test.protein.score < best_score:
                best_protein = copy.deepcopy(greedy_test.protein)
                best_score = best_protein.score
    
    # stop the timer and append time of the split, with the split number and the best score to total lists
    end_time = time.time()
    total_time.append(end_time - start_time)
    splits_total.append(key)
    best_scores.append(best_score)

print(best_protein)
print(best_score)

# create a dataframe with split numbers, scores and befores for every iteration, for a distribution of scores
df_exp_greedy_big_parts['split_numbers'] = splits
df_exp_greedy_big_parts['scores'] = scores
df_exp_greedy_big_parts['befores'] = befores

# create a dataframe with splits, best scores and total time for the summarized results of one experiment
df_exp_greedy_big_parts_total['split_number'] = splits_total
df_exp_greedy_big_parts_total['best_scores'] = best_scores
df_exp_greedy_big_parts_total['total_time'] = total_time

# visualize the best protein
visualize_protein(best_protein, 3)

print(df_exp_greedy_big_parts.head())
#df_exp_greedy_big_parts.to_csv(path_or_buf=r'C:\Users\alaya\Documents\Data_Greedy\df_exp_greedy_big_parts')

print(df_exp_greedy_big_parts_total.head())
#df_exp_greedy_big_parts_total.to_csv(path_or_buf=r'C:\Users\alaya\Documents\Data_Greedy\df_exp_greedy_big_parts_total')