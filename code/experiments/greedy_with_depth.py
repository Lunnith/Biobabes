from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms import greedy as gr
import copy
import pandas as pd
import time    

# select a sequence and create two empty dataframes for all scores and for summary results
#sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
sequence = "HHPHHHPHPHHHPH"
df_exp_greedy = pd.DataFrame()
df_exp_greedy_total = pd.DataFrame()

# make empty lists to fill during the experiment
splits = []
splits_total = []
scores = []
best_scores = []
befores = []
total_time = []

# make a dictionary with the splits as keys and the number of iterations per 'before' as the value
dict_splits = {1: 10, 2: 5, 3: 4, 4: 3}
#{1: 2000, 2: 1000, 3: 667, 4: 500, 5: 400} # this one takes about 3 hours
for key, value in dict_splits.items():
    
    # start the clock for measuring the time of one every split
    start_time = time.time()

    # set best score to 1
    best_score = 1

    # print the key to see where you are in the experiment
    print(key)

    # for every before in the range of the split do the number of iterations for that split
    for before in range(key):

        for i in range(value):

            # make a protein, initialize greedy with the right parameters and run the algorithm
            protein = Protein(sequence)
            greedy_test = gr.Greedy(protein, 2, splits = key, before = before)
            greedy_test.run()

            # append for each iterations the score, the before and the split
            scores.append(greedy_test.protein.score)
            befores.append(before)
            splits.append(key)

            # if the score is lower than the best score, set best score at protein score
            if greedy_test.protein.score < best_score:
                best_protein = copy.deepcopy(greedy_test.protein)
                best_score = greedy_test.protein.score
            print(greedy_test.states)
    # stop the timer and append time of the split, with the split number and the best score to total lists
    end_time = time.time()
    total_time.append(end_time - start_time)
    splits_total.append(key)
    best_scores.append(best_score)
    if best_score <= best_scores[-1]:
        best_protein_total = best_protein

# create a dataframe with split numbers, scores and befores for every iteration, for a distribution of scores
df_exp_greedy['split_numbers'] = splits
df_exp_greedy['scores'] = scores
df_exp_greedy['befores'] = befores

# create a dataframe with splits, best scores and total time for the summarized results of one experiment
df_exp_greedy_total['split_number'] = splits_total
df_exp_greedy_total['best_scores'] = best_scores
df_exp_greedy_total['total_time'] = total_time

# check the data and visualize the best protein
print(df_exp_greedy.head())
print(df_exp_greedy_total.head())
visualize_protein(best_protein_total, 2)

# plot the data using a boxplot for ever split to show score distributions
sns.boxplot(data = df_exp_greedy, x = 'split_numbers', y = 'scores')
ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('Split')
plt.ylabel('Score')
plt.title('Score distributions Greedy with depth 2000 iterations')
plt.show()