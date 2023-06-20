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

    plt.title("Hill_climber")
    plt.legend(range(1, max_n))
    plt.xlim(left=0, right=len(sequence))
    # plt.ylim()
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

    plt.title("Randomise", fontweight='bold')
    plt.xlabel("Length protein-sequence", loc='right')
    plt.ylabel("Time per 50 iterations\n(in seconds)", loc='top')
    plt.show()


# algorithm_timing(sequence)

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
        
        important_parts = ImportantParts(protein, 3)
        important_parts.run_in_parts(n=1000, split_on_P=True, split_on_size=False)

        end = time.time()
        times[seq_length] = (end-start)

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
    plt.show()


algorithm_timing(sequence)

########################### Greedy
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
        
        greedy_test = gr.Greedy(protein, 3, 2)
        greedy_test.run_k()

        end = time.time()
        times[seq_length] = (end-start)

        #Make safety net for exponential algorithms
        if end-start > 1200: #That's 20 minutes
            break

    total_end = time.time()
    print("Total runtime:", total_end-total_start)
    plt.plot(times.keys(), times.values())
    plt.xlim(left=0, right=len(sequence))

    plt.title("Greedy", fontweight='bold')
    plt.xlabel("Length protein-sequence", loc='right')
    plt.ylabel("Time per run\n(in seconds)", loc='top')
    plt.show()


# algorithm_timing(sequence)