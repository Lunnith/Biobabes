from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
from ..algorithms.hill_climber import Hill_climber
from ..algorithms.simulated_annealing import SimulatedAnnealing
from ..visualization.visualize import *
from ..algorithms.randomise import random_assignment
import copy
import time


def unpack_scores(scores, improvement):
    all_scores = {}
    for iteration in range(iterations+1):
        all_scores[iteration] = []
        if improvement[iteration] == "Y" or improvement[iteration] == "S":
            improved_score = scores[iteration]
            all_scores[iteration].append(improved_score)
        else:
            all_scores[iteration].append(None)
    return all_scores

for test in range(0,5):
    sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    protein = random_assignment(Protein(sequence), 3)
    while len(sequence) != len(protein.sequence_list) or protein.score < -10:
        protein = random_assignment(Protein(sequence), 3)
    testing_protein = copy.deepcopy(protein)

    # For these experiments, we'll keep the starting n at 10, as this n decreases during the run
    # The best value for n will be found in the hill_climber vs sim_annealing experiment, 
    # which will be run with the best parameters found in this experiment.

    # Define repeating testing params:
    iterations = 2500
    dimensions = 3
    folded = True
    start_n = 10

    # define variable testing params:
    testing_temperatures = range(1, 25+1)

    #Try different starting temperatures and different temperature schemes:
    start = time.time()

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    ax = 0
    for scheme in ['exp', 'lin']:
        results_dict = {}

        for temp in testing_temperatures:
            sa = SimulatedAnnealing(protein, start_n, folded=folded, dimensions=dimensions, temperature=temp, temp_scheme=scheme)
            sa_protein, sa_lowest_score, sa_scores, sa_improvement = sa.run_i_iterations(testing_protein, iterations, start_n, sim_annealing=True, sample_number=temp)
            sa_scores = unpack_scores(sa_scores, sa_improvement)

            temp_df = pd.DataFrame.from_dict(sa_scores, orient='index')
            temp_df = temp_df.fillna(method='ffill')
            temp_df['Average'] = temp_df.mean(axis=1)

            results_dict[temp] = temp_df["Average"]
        results = pd.DataFrame.from_dict(results_dict, orient='columns')

        sns.lineplot(data=results, hue_order=results.columns, dashes=False, ax=axs[ax])
        ax += 1

    for i in range(ax):
        axs[i].set_xlabel("Iterations")
        axs[i].set_ylabel("Score")
        axs[i].set_xlim(0, iterations)
        axs[i].legend(title="Starting temperature", ncol=5, loc='upper center')
        axs[i].set_title(f"Temperature scheme = {['Exponential', 'Linear'][i]}")

    end = time.time() 
    print(f"runtime {end-start} seconds")   

    plt.savefig(f"{iterations} Iterations {test}")
    plt.clf()
