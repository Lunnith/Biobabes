from code.classes.protein import Protein

from code.algorithms.randomise import *
from code.algorithms.greedy import Greedy
from code.algorithms.depth_first import DepthFirst
from code.algorithms.important_parts import ImportantParts
from code.algorithms.hill_climber import Hill_climber
from code.algorithms.simulated_annealing import SimulatedAnnealing

import time
import pandas as pd

df_short_scores = pd.DataFrame()
df_long_scores = pd.DataFrame()

test_protein_short = Protein("HHPHHHPHPHHHPH")
test_protein_long = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

# ---------------------------- Random ----------------------------
random_protein_short, score_list = random_reassignment(test_protein_short, 3, k=100000)
random_protein_long = random_reassignment(test_protein_long, 3, k=100000)

df_short_scores["Random algorithm"] = random_protein_short.score
df_long_scores["Random algorithm"] = random_protein_long.score

# ---------------------------- Greedy ----------------------------
greedy_scores_short = []
greedy_scores_long = []

for i in range(25000):
    test_protein_short = Protein("HHPHHHPHPHHHPH")

    greedy_protein_short = Greedy(test_protein_short, 3, splits=1)
    greedy_protein_short.run()
    greedy_scores_short.append(greedy_protein_short.protein.score)

for i in range(25000):
    test_protein_long = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

    greedy_protein_long = Greedy(test_protein_long, 3, splits=1)
    greedy_protein_long.run()
    greedy_scores_long.append(greedy_protein_long.protein.score)

df_short_scores["Greedy algorithm"] = min(greedy_scores_short)
df_long_scores["Greedy algorithm"] = min(greedy_scores_long)

# ---------------------------- Greedy with beam search ----------------------------
greedy_beam_scores_short = []
greedy_beam_scores_long = []

for i in range(500):
    for before in range(4):
        test_protein_short = Protein("HHPHHHPHPHHHPH")

        greedy_beam_protein_short = Greedy(test_protein_short, 3, splits=4, before=before)
        greedy_beam_protein_short.run()
        greedy_beam_scores_short.append(greedy_beam_protein_short.protein.score)

for i in range(500):
    for before in range(4):
        test_protein_long = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

        greedy_beam_protein_long = Greedy(test_protein_long, 3, splits=4, before=before)
        greedy_beam_protein_long.run()
        greedy_beam_scores_long.append(greedy_beam_protein_long.protein.score)

df_short_scores["Greedy algorithm with depth"] = min(greedy_beam_scores_short)
df_long_scores["Greedy algorithm with depth"] = min(greedy__beam_scores_long)

# ---------------------------- Depth First with P and directions pruning ----------------------------
depth_first_protein = DepthFirst(test_protein_short, 3)
depth_first_protein.run(P_pruning=True, directions_pruning=True)

df_short_scores["Depth First algorithm"] = depth_first_protein.protein.score

# ---------------------------- Depth First - Important Parts ----------------------------
# HIER BEZIG --> MAAK AF EN VOEG OOK TIJD TOE
important_parts_scores_short_on_P = []
important_parts_scores_long_on_P = []

important_parts_scores_short_on_size = []
important_parts_scores_long_on_size = []

important_parts_protein_short_on_P = ImportantParts(test_protein_short, 3)
important_parts_protein_short_on__P.run(iterations=1000, split_on_P=True, split_on_size=False)

important_parts_protein_short_on_P = ImportantParts(test_protein_short, 3)
important_parts_protein_short_on__P.run(iterations=1000, split_on_P=True, split_on_size=False)



# ---------------------------- Hill Climber ----------------------------
hill_climber_protein = Hill_climber(test_protein_A, 3)
hill_climber_protein.run_i_iterations(test_protein_A, iterations=1000, bonds=1)

print(f'Value of the folding after Hill Climber:'
      f'{hill_climber_protein.protein.score}')

# ---------------------------- Simulated Annealing ----------------------------
simanneal_protein = SimulatedAnnealing(test_protein_I, 10, False, 3,  temperature=20)
simanneal_protein.run_i_iterations(test_protein_I, iterations=1000, bonds=10)

print(f'Value of the folding after Simulated Annealing:'
        f'{simanneal_protein.protein.score}')