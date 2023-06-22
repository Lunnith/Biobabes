from ..classes.protein import Protein
from .randomise import random_assignment
import matplotlib.pyplot as plt
import pandas as pd
import random
import time
import copy

class Hill_climber():
    """
    The Hill_climber algorithm is meant to take a protein and make small changes to hopefully improve 
    the score of this protein. If the score of the protein did improve, it keeps the adjustment.

    By repeating this process over and over again, the Hill-climber finds an optimum score.
    However, by keeping each change that improved the score, it could be possible that an even better 
    improvement is looked over. This results in the Hill climber getting stuck in a local optimum.
    ...
    Attributes:
    ----------
    protein: Protein
        The protein to improve.
    dimensions: int
        Number of dimensions in which the protein can fold.
    prints: bool
        Wether or not to print small updates.
    folded: bool
        Wether or not the given protein has already been folded.
    lowest_score: int
        Best score found yet.
    directions: dict
        Possible directions, depending on the amount of dimensions used

    Methods:
    change_bond(protein, given_index=None, skip_bonds=None):
        Change one single bond
    refold(protein, index):
        Refold protein from updated bond
    check_validity(protein):
        Check wether the protein has folded over itself
    refold_into_valid_state(protein):
        Refold protein into a valid state
    change_n_bonds(protein, n):
        Change n different bonds
    check_score(protein):
        Check the score of a given protein
    check_solution(new_protein):
        Check if the new protein has a better score
    run_i_iterations(protein, iterations, bonds):
        Run an amount of iterations of changes
    plot_hillclimb(iterations, scores, n):
        Plot the climb made
    optimize_graph(n, iterations):
        Give plot a nice appearance and show the plot
    experiment(protein, iterations, sample_size=1, max_n=10):
        Run an experiment with multiple runs per different amount of bonds to change
    """
    def __init__(self, protein: Protein, dimensions=3, prints=False, folded=False) -> None:
        #Initiate first folding
        if folded: self.protein = protein
        else: 
            if len(protein.sequence_list) > 1:
                raise Exception("By default, Hill_climber assumes that you input an unfolded protein.\n\
                If you want to input a protein that is already folded, use the keyword 'folded=True'")
            self.protein = random_assignment(protein, dimensions)
        self.lowest_score = self.protein.score

        if prints == True:
            self.prints = True
        else: self.prints = False
        
        #Set direction options for dimensions
        if dimensions == 2:
            self.directions = {1: [1, 0, 0, 1], -1: [-1, 0, 0, -1], 2: [0, 1, 0, 2], -2: [0, -1, 0, -1]}
        if dimensions == 3:
            self.directions = {1: [1, 0, 0, 1], -1: [-1, 0, 0, -1], 2: [0, 1, 0, 2], -2: [0, -1, 0, -2], 3: [0, 0, 1, 3], -3: [0, 0, -1, -3]}


    def change_bond(self, protein: Protein, given_index=None, skip_bonds=None) -> tuple[Protein, int]:
        """
        This function takes a folded protein and optionally the index of the bond to change.
        If no index is given, a random bond of the protein is changed.
        The change in direction is randomly chosen.

        Note:
        This function only changes the direction of the bond and doesn't update any coördinates!
        """
        if given_index == None:
            index_changing_bond = random.randint(2, len(protein.sequence_list)-1)
            if skip_bonds != None: #If there are bonds that we want to ignore
                while index_changing_bond in skip_bonds:            
                    index_changing_bond = random.randint(2, len(protein.sequence_list)-1)
        else:
            index_changing_bond = given_index

        acid = protein.sequence_list[index_changing_bond]
        acid.location_valid = False

        tried_directions = set()
        tried_directions.add(acid.step)

        while acid.location_valid == False:
            new_step = random.sample(self.directions.keys(), 1)[0]
            tried_directions.add(new_step)

            #Overwrite step of previous acid
            protein.sequence_list[index_changing_bond - 1].step = new_step
            #fix format of directions
            direction = tuple(self.directions[new_step])

            #Create new bond based on new direction
            protein.create_bond(acid, protein.sequence_list[index_changing_bond - 1], direction)

            # if protein can't fold anymore, return shorter folded protein and a False index
            if 0 in tried_directions: #Safety net for last aminoacid
                tried_directions.remove(0)
            if tried_directions == self.directions.keys():
                if self.prints: print('change_bond has gotten stuck:', end=" ")
                return protein, False     
            
        return protein, index_changing_bond
    

    def refold(self, protein: Protein, index: int) -> Protein:
        """
        This function takes a protein where a bond has been redirectioned 
        and also takes the index of this bond. Then, it updates the coördinates of the aminoacids.

        Note:
        This function can make the protein fold over itself and leaves it without hh/ch/cc bonds!
        """
        #reset the old hh/ch/cc bonds
        protein.hh_ch_bonds = []
        protein.cc_bonds = []

        #Retrieve folded protein up until the changed bond
        used_coords = set()
        for acid in range(0, index):
            used_coords.add(tuple(protein.sequence_list[acid].location))
            protein.used_coordinates = used_coords

        #Refold from the changed bond
        for acid in range(index, len(protein.sequence_list)):
            direction = self.directions[protein.sequence_list[acid-1].step]
            protein.create_bond(protein.sequence_list[acid], protein.sequence_list[acid-1], direction)

        return protein
         

    def check_validity(self, protein: Protein) -> list:
        """
        This function takes a (re)folded protein and checks wether or not it has folded over itself.
        Then, it returns a list of coördinates where multiple aminoacids appear
        """
        self.double_coords = []
        used_coords = set()

        for acid in protein.sequence_list:
            if tuple(acid.location) not in used_coords: #Keep track of all locations in use
                used_coords.add(tuple(acid.location))
            else: #If already in use, add to list of doubles
                index = protein.sequence_list.index(acid)
                self.double_coords.append(index)

        return self.double_coords
    

    def refold_into_valid_state(self, protein: Protein) -> Protein:
        """
        This function refolds the protein, 
        checks it it has folded over itself and fixes this if nessecary.

        If this was possible, it returns the new protein.
        If the protein could not fold into a valid state, it returns False
        """        
        #Fix the aminoacids that have been folded over eachother
        while len(self.check_validity(protein)) > 0:

            # Do this in chronological order
            protein, changed_bond = self.change_bond(protein, given_index=self.double_coords[0])
            if changed_bond == False: #If protein could not fold into a valid state
                if self.prints: print(f"The changes resulted in an unfixable folding version")
                return False
            
            self.refold(protein, changed_bond)
        return protein
    

    def change_n_bonds(self, protein: Protein, n: int="self_n") -> Protein:
        """
        Loop n times over the function change_bond,
        so that n bonds have been changed.
        Then, refold the protein into a valid state
        """
        if n != "self_n":
            self.n = n
        for i in range(self.n):
            changed_bonds = []
            continued = 0

            protein, changed_bond = self.change_bond(protein, skip_bonds=changed_bonds) #you don't want to change the same bond multiple times
            if changed_bond == False: #If protein could not fold into a valid state, skip this change
                continued += 1
                if self.prints: print(f"Changing {continued} less bonds")
                continue
            else: 
                protein = self.refold(protein, changed_bond)
                changed_bonds.append(changed_bond)

        if continued == self.n: #If all n changes were invalid
            if self.prints: print(f"All {self.n} initiated changes in bonds were invalid.")
            return False

        protein = self.refold_into_valid_state(protein)
        if protein == False: #If protein could not refold into a valid state
            return False

        return protein


    def check_score(self, protein: Protein) -> None:
        """
        This function checks the score of a given protein.
        """
        protein.score = 0
        for acid in range(len(protein.sequence_list)):
            protein.sequence_list[acid].check_interactions(protein, index=acid+1)
    

    def check_solution(self, new_protein: Protein) -> None:
        """
        Checks if the input protein has a better score than the best protein yet.
        """
        if new_protein.score < self.lowest_score:
            self.improvement.append("Y") #Yes
            self.lowest_score = new_protein.score
            self.protein = new_protein
            if self.prints:
                print(f"Score updated to", new_protein.score, end=" ")
                printed = True
                return printed

        elif new_protein.score == self.lowest_score:
            self.improvement.append("S") #Same
        else:
            self.improvement.append("N") #No


    def run_i_iterations(self, protein: Protein, iterations: int, bonds: int) -> tuple[Protein, int, list, list]:
        """
        runs the change_n_bonds for a given amount of iterations.
        Then, returns the protein with the best score, the actual score, the list of scores that 
        were considered and a list with on the index of each score, wether or not it was an improvement.
        """
        self.iterations = iterations
        starting_score = protein.score
        scores = [starting_score]
        self.improvement = ["Y"]
        self.n = bonds

        self.lowest_score = starting_score
        for i in range(iterations):
            protein = copy.deepcopy(self.protein)
            new_protein = self.change_n_bonds(protein)
            if new_protein == False: #If change turned out invalid, skip this change
                if self.prints: (f"Skipping iteration {i}, change turned out invalid")
                scores.append(None)
                self.improvement.append("NaN")
                continue

            self.check_score(new_protein)
            scores.append(new_protein.score)
            updated = self.check_solution(new_protein)
            if updated == True: print(f"in iteration {i}.")

        return self.protein, self.lowest_score, scores, self.improvement
    

    def plot_hillclimb(self, iterations: int, scores: list, n: int) -> None:
        """
        Plots only the improvement graph for this specific n.
        If a plt object is already existing, this function plots over it.
        Does not do any appearance changes to the plot, nor shows it.
        """
        colors = ['firebrick', 'orangered', 'darkorange','gold', 'yellowgreen', 'lightgreen', 'turquoise', 'lightskyblue', 'plum', 'lightpink']
        plt.plot(iterations, scores, "-", linewidth=2, c=colors[n-1], label=n)


    def optimize_graph(self, n: int, iterations: int) -> None:
        """
        Uses the existing plt object where all hill-climbs have been plotted,
        Makes it aesthetically pleasing and then shows it.
        """
        plt.xlim(left=0, right=iterations)
        plt.ylim(top=0)
        plt.title("Hill-climber", fontweight='bold')
        plt.xlabel("Iterations", loc='right')
        plt.ylabel("Score", loc='top')
        plt.legend(range(1, n+1), title="Minimum amount of bonds changed per iteration", ncol=n//2)
        plt.show()


    def experiment(self, protein: Protein, iterations: int, sample_size=1, max_n=10) -> Protein:
        """
        Runs an experiment with a given protein. The sample size is the amount of times to run the algorithm.
        It then runs the algorithm for each n amount of bonds to change, with the given amount of iterations.
        Afterwards, it takes the average improvement made in all samples of runs,
        and plots the improvement for each n per iteration.
        """
        if self.prints: print("\n\nStarting score =", protein.score, "Length protein =", len(self.protein.sequence_list))
        
        start = time.time()
        if max_n > 10: raise ValueError("Please only insert a max_n of 10 or smaller")

        first_random_fold = copy.deepcopy(protein)
        best_protein = first_random_fold
        best_score = first_random_fold.score
        
        for n in range(1, max_n+1):
            if self.prints: print(f"\nStarting new N: {n}")

            all_scores = {}
            for iteration in range(iterations+1):
                all_scores[iteration] = []

            for sample_run in range(sample_size):
                protein_for_n, lowest_score_for_n, scores, improvement = self.run_i_iterations(first_random_fold, iterations, n)

                if lowest_score_for_n < best_score:
                    best_protein = protein_for_n
                    best_score = lowest_score_for_n

                for iteration in range(iterations+1):
                    if improvement[iteration] == "Y" or improvement[iteration] == "S":
                        improved_score = scores[iteration]
                        all_scores[iteration].append(improved_score)
                    else:
                        all_scores[iteration].append(None)
                        
            temp_df = pd.DataFrame.from_dict(all_scores, orient='index')
            temp_df = temp_df.fillna(method='ffill')
            temp_df['Average'] = temp_df.mean(axis=1)
            self.plot_hillclimb(temp_df.index, temp_df['Average'], n)
        
        end = time.time()
        if self.prints: print(f"\n\nBest score has become", best_score)
        if self.prints: print(f"Runtime experiment: {end-start} seconds.\n")
        self.optimize_graph(max_n, iterations)

        return best_protein