from ..algorithms.randomise import random_reassignment
import matplotlib.pyplot as plt
from ..visualization.visualize import *
from ..classes.protein import Protein

import time

# start the timer
start_time = time.time()

# set the sequence and initiate the protein
sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
protein = Protein(sequence)

# call random reassingment with a number of iterations andd the protein and return the best protein, all the scores and the best score
k = 10000
best_protein, scores, best_score = random_reassignment(protein, 3, k)

# end the timer and determine the total time
end_time = time.time()
total_time = round((end_time - start_time), 2)

# print and visualize summerizing results
print('time', total_time)
print('best score', best_score)
visualize_protein(best_protein, 3)

# visualize a histogram with scores and also show the best score and the used time
plt.hist(scores, rwidth = 0.75, bins = 22, align = 'right') # adjust the bins depending on the results to get the best visualization
plt.title(f'Random algoritme {k} iteraties')
plt.xlabel('Score')
plt.ylabel('Frequentie')
plt.text(x = -20, y = 1200, s = f'Best score: {best_score}, time: {total_time} sec', fontsize = 10) # adjust placement of the text depending on results
plt.show()
