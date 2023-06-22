from code.classes.protein import Protein
from code.visualization.visualize import *
from code.algorithms.hill_climber import *
from code.algorithms.randomise import *
from code.algorithms.important_parts import *
from code.algorithms import greedy as gr
from operator import add
import matplotlib.pyplot as plt
import random
import time

sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
# sequence = "HCPHPCPHPCHCH"


######################## Hill_climber
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

    plt.xlim(left=0, right=len(sequence))
    plt.title("Hill-climber", fontweight='bold')
    plt.xlabel("Length protein-sequence", loc='right')
    plt.ylabel("Time per 500 iterations\n(in seconds)", loc='top')
    plt.legend(range(1, n+1), title="Minimum amount of bonds changed per iteration", ncol=n//2)
    plt.show()

# timing_for_multiple_n(sequence)

################################ Randomise
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
        
        start = time.time()
        folded_protein, score_list = random_reassignment(protein, 3, k=50)
        end = time.time()
        times[seq_length] = (end-start)

        #Make safety net for exponential algorithms
        if end-start > 1200: #That's 20 minutes
            break

    total_end = time.time()
    print("Total runtime:", total_end-total_start)
    plt.plot(times.keys(), times.values())
    plt.xlim(left=0, right=len(sequence))


def timing_for_multiple_iteration_amounts(sequence, max_n=500):
    """
    
    """
    start = time.time()
    for n in range(50, max_n+50, 50):
        algorithm_timing(sequence)
    end = time.time()
    print("Runtime multiple n:", end-start)

    plt.xlim(left=0, right=len(sequence))
    plt.title("Randomise", fontweight='bold')
    plt.xlabel("Length protein-sequence", loc='right')
    plt.ylabel("Time per n iterations\n(in seconds)", loc='top')
    plt.legend(range(50, max_n+50, 50), title="Amount of iterations", ncol=(n/50)//2)
    plt.show()

# timing_for_multiple_iteration_amounts(sequence)

########################## Depth First

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
        
        start = time.time()
        
        depth_first = DepthFirst(protein, 3)
        depth_first.run(P_pruning=False, directions_pruning=False)

        end = time.time()
        times[seq_length] = (end-start)
        print(seq_length, "Length done")

        #Make safety net for exponential algorithms
        if end-start > 1200: #That's 20 minutes
            break

    total_end = time.time()
    print("Total runtime:", total_end-total_start)
    plt.plot(times.keys(), times.values())
    plt.xlim(left=0, right=len(sequence))

    plt.title("Depth First", fontweight='bold')
    plt.xlabel("Length protein-sequence", loc='right')
    plt.ylabel("Time per run\n(in seconds)", loc='top')

    plt.savefig("DepthFirst_test")
    plt.show()


algorithm_timing(sequence)

########################### Greedy
def algorithm_timing(sequence, iterations, split):
    """

    """
    total_start = time.time()
    times = {}
    for seq_length in range(len(sequence)):
        times[seq_length] = None

    for seq_length in range(len(sequence)):
        if seq_length <2:
            continue
        
        
        start = time.time()
        for i in range(iterations):
            protein = Protein(sequence[:seq_length+1])
            greedy_test = gr.Greedy(protein, 3, split)
            greedy_test.run_k()

        end = time.time()
        times[seq_length] = (end-start)

        #Make safety net for exponential algorithms
        if end-start > 1200: #That's 20 minutes
            break

    total_end = time.time()
    print("Total runtime:", total_end-total_start)
    plt.plot(times.keys(), times.values())
    plt.savefig(f"Greedy-test_run2_{iterations}")


def timing_for_multiple_splits(sequence, iterations=20):
    """
    
    """
    start = time.time()
    splits = [1, 2, 3, 4, 5]
    for split in splits:
        algorithm_timing(sequence, iterations, split)
    end = time.time()
    print("Runtime multiple n:", end-start)

    plt.title("Greedy", fontweight='bold')
    plt.xlim(left=0, right=len(sequence))
    plt.xlabel("Length protein-sequence", loc='right')
    plt.ylabel(f"Time per {iterations} iterations\n(in seconds)", loc='top')
    plt.legend(range(1, len(splits)+1), title="Splits", ncol=len(splits)//2)
    # plt.show()
    plt.savefig(f"Greedy-test_run2_{iterations}")
    plt.clf()

# timing_for_multiple_splits(sequence, iterations=10)
# timing_for_multiple_splits(sequence)
# timing_for_multiple_splits(sequence, iterations=30)

# timing_for_multiple_splits(sequence, iterations=40)
# timing_for_multiple_splits(sequence, iterations=50)
