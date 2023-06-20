from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.hill_climber import *
import time

#NOTE TO SELF: 2D version is buggy
# protein = Protein("HPCPHHPCCHPHCCHPHHHCCCPPHCPHCPHCCCPHHHCPPCCHPCCCHPHCCHHHHPCCCPPPCH")
# protein = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
# sequence = "HCPHCH"

#Creating time thingy
# def algorithm_timing(sequence, algorithm, *kwargs):
def algorithm_timing(sequence):
    """
    
    """
    total_start = time.time()
    times = {}
    for seq_length in range(len(sequence)):
        times[seq_length] = None

    for seq_length in range(len(sequence)):
        if seq_length <2:
            continue
        protein = Protein(sequence[:seq_length+1])
        hill_climber = Hill_climber(protein)

        start = time.time()
        hill_climber.run_n_iterations(protein, 500, 1)
        end = time.time()
        times[seq_length] = (end-start)
        
        #Make safety net for exponential algorithms
        if end-start > 1200: #That's 20 minutes
            break

    total_end = time.time()
    print("Total runtime:", total_end-total_start)
    plt.plot(times.keys(), times.values())
    # plt.title("Hill_climber")
    # plt.show()

# algorithm_timing(sequence)

def timing_for_multiple_n(sequence, max_n=10):
    """
    
    """
    start = time.time()
    for n in range(1, max_n):
        algorithm_timing(sequence)
    end = time.time()
    print("Runtime multiple n:", end-start)

    plt.title("Hill_climber")
    plt.legend(range(1, max_n))
    plt.xlim(left=0, right=len(sequence))
    # plt.ylim()
    plt.title("Hill-climber", fontweight='bold')
    plt.xlabel("Length protein-sequence", loc='right')
    plt.ylabel("Time per 500 iterations\n(in seconds)", loc='top')
    plt.legend(range(1, n+1), title="Minimum amount of bonds changed per iteration", ncol=n//2)
    plt.show()

timing_for_multiple_n(sequence)

# # protein = Protein("HPPCCHPPH")

# print("Starting score =", hill_climber.lowest_score)
# # hill_climber.change_one_bond(protein)
# # visualize_protein(folded_protein, 3)



# folded_protein, best_score = hill_climber.run_n_iterations(protein, 1000, 1)
# print(best_score)
# visualize_protein(folded_protein, 3)
# if folded_protein.score > -20:
#     best_protein = hill_climber.experiment(folded_protein, 500, 50, max_n=10)
#     visualize_protein(best_protein, 3)
# else: print("Score too good:", folded_protein.score)

