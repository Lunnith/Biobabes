from ..algorithms.randomise import random_reassignment
import matplotlib.pyplot as plt
from ..visualization.visualize import *
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
import time

# # create an empty dataframe
# df_random = pd.DataFrame()

# # start the timer
# start_time = time.time()

# # set the sequence and initiate the protein
# sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
# protein = Protein(sequence)

# # call random reassingment with a number of iterations andd the protein and return the best protein, all the scores and the best score
# best_protein, scores, best_score = random_reassignment(protein, 3, k = 100000)

# # end the timer and determine the total time
# end_time = time.time()
# total_time = (end_time - start_time)

# # print and visualize interesting results
# print(total_time)
# print('best score', best_score)
# visualize_protein(best_protein, 3)

# df_random['scores'] = scores

# df_random.to_csv(path_or_buf=r'C:\Users\alaya\Documents\Data_Greedy\df_exp_random_100000i')

#data_random_100000i = pd.read_csv(r'C:\Users\alaya\Documents\Data_Greedy\df_exp_random_100000i.csv')
data_random_10000i = pd.read_csv(r'C:\Users\alaya\Documents\Data_Greedy\df_exp_random_10000i.csv')
# print(len(data_random_20mil))
# print(data_random_20mil.head())
scores = data_random_10000i['scores']
print(scores.min())
# plt.hist(data_random_100000i['scores'], rwidth = 0.75, bins = 26)
# plt.show()

plt.hist(data_random_10000i['scores'], rwidth = 0.75, bins = 22, align = 'right')
plt.title('Random algoritme 10.000 iteraties')
plt.xlabel('Score')
plt.ylabel('Frequentie')
plt.text(x = -18, y = 1200, s = 'Beste score: -22', fontsize = 10)
plt.show()
