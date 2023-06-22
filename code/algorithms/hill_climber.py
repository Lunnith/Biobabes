from .randomise import random_assignment
import random
import time
import copy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

class Hill_climber():
    def __init__(self, protein, dimensions=3, prints=False):
        self.protein = protein
        self.double_coords = []
        if prints == True:
            self.prints = True
        else: self.prints = False
        
        #Set options for dimensions
        if dimensions == 2:
            self.directions = {1: [1, 0, 0, 1], -1: [-1, 0, 0, -1], 2: [0, 1, 0, 2], -2: [0, -1, 0, -1]}
        if dimensions == 3:
            self.directions = {1: [1, 0, 0, 1], -1: [-1, 0, 0, -1], 2: [0, 1, 0, 2], -2: [0, -1, 0, -2], 3: [0, 0, 1, 3], -3: [0, 0, -1, -3]}

        #Initiate first random folding
        self.protein = random_assignment(protein, dimensions)
        self.lowest_score = protein.score

    
    def change_bond(self, protein, given_index=None, skip_bonds=None):
        """
        This function takes a folded protein and optionally the index of the bond to change.
        If no index is given, a random bond of the protein is changed.
        The change in direction is randomly chosen.

        Note:
        This function only changes the direction of the bond and doesn't update any coördinates!
        """
        if given_index == None:
            index_changing_bond = random.randint(2, len(protein.sequence_list)-1)
            if skip_bonds != None:
                while index_changing_bond in skip_bonds:            
                    index_changing_bond = random.randint(2, len(protein.sequence_list)-1)
        else:
            index_changing_bond = given_index

        old_step = protein.sequence_list[index_changing_bond].step
        new_step = old_step
        acid = protein.sequence_list[index_changing_bond]

        #Determine new step
        tried_directions = set()
        tried_directions.add(old_step)
        acid.location_valid = False
        
        while acid.location_valid == False:
            new_step = random.sample(self.directions.keys(), 1)[0]
            tried_directions.add(new_step)

            #Overwrite step of previous acid
            protein.sequence_list[index_changing_bond - 1].step = new_step
            #fix format of directions
            direction = tuple(self.directions[new_step])

            #Create new bond based on new direction
            protein.create_bond(acid, protein.sequence_list[index_changing_bond - 1], direction)

            # if protein can't fold anymore, return shorter folded protein
            if 0 in tried_directions: #Safety net for last aminoacid
                tried_directions.remove(0)
            if tried_directions == self.directions.keys():
                if self.prints: print('change_bond has gotten stuck:', end=" ")
                return protein, False     
            
        return protein, index_changing_bond
    
    def refold(self, protein, index):
        """
        This function takes a protein where a bond has been redirectioned 
        and also takes the index of this bond. Then, it updates the coördinates of the aminoacids.

        Note:
        This function can make the protein fold over itself!
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
         

    def check_score(self, protein):
        """
        This function checks the score of a given protein
        """
        protein.score = 0
        for acid in range(len(protein.sequence_list)):
            protein.sequence_list[acid].check_interactions(protein, index=acid+1)


    def check_validity(self, protein):
        """
        This function takes a (re)folded protein and checks wether or not it has folded over itself.
        Then, it returns a list of coördinates where multiple aminoacids appear
        """
        self.double_coords = []
        used_coords = set()
        for acid in protein.sequence_list:
            if tuple(acid.location) not in used_coords:
                used_coords.add(tuple(acid.location))
            else:
                index = protein.sequence_list.index(acid)
                self.double_coords.append(index)

        return self.double_coords
    

    def refold_into_valid_state(self, protein):
        """
        This function refolds the protein, 
        checks it it has folded over itself and fixes this if nessecary.

        If this was possible, it returns the new protein.
        If the protein could not fold into a valid state, it returns False
        """        
        #Fix the aminoacids that have been folded incorrectly
        while len(self.check_validity(protein)) > 0:

            # Fix this in chronological order
            protein, changed_bond = self.change_bond(protein, given_index=self.double_coords[0])
            if changed_bond == False: #If protein could not fold into a valid state
                if self.prints: print(f"The changes resulted in an unfixable folding version")
                return False
            
            self.refold(protein, changed_bond)
        return protein


    def change_n_bonds(self, protein, n):
        """
        Loop n times over the function change_bond,
        so that n bonds have been changed.
        Then, refold the protein into a valid state
        """
        for i in range(n):
            continued = 0
            changed_bonds = []

            protein, changed_bond = self.change_bond(protein, skip_bonds=changed_bonds) #you don't change the same bond multiple times
            if changed_bond == False: #If protein could not fold into a valid state, skip this try
                continued += 1
                if self.prints: print(f"Changing {continued} less bonds")
                continue
            else: 
                protein = self.refold(protein, changed_bond)
                changed_bonds.append(changed_bond)

        if continued == n: #If all n changes were invalid
            if self.prints: print(f"All {n} changes were invalid.")
            return False

        protein = self.refold_into_valid_state(protein)
        if protein == False: #If protein could not refold into a valid state
            return False

        return protein



    def run_n_iterations(self, protein, iterations, bonds):
        """
        runs the change_n_bonds
        """
        starting_score = protein.score
        self.lowest_score = starting_score
        scores = [starting_score]
        improvement = ["Y"]

        for n in range(iterations):
            protein = copy.deepcopy(self.protein)
            new_protein = self.change_n_bonds(protein, bonds)
            if new_protein == False: #If change turned out invalid, skip this change
                if self.prints: (f"Skipping iteration {n}, change turned out invalid")
                scores.append(None)
                improvement.append("NaN") #NaN
                continue

            self.check_score(new_protein)
            scores.append(new_protein.score)


            ###### Make this its own function
            if new_protein.score < self.lowest_score:
                if self.prints: print(f"Iteration {n}: Score updated to", new_protein.score)
                improvement.append("Y") #Yes
                self.lowest_score = new_protein.score
                self.protein = new_protein

            elif new_protein.score == self.lowest_score:
                improvement.append("S") #Same
            else:
                improvement.append("N") #No

        return self.protein, self.lowest_score, scores, improvement
    

    def plot_hillclimb(self, iterations, scores, n):
        """
        Roughly visualises the improvement of the algorithm
        """
        colors = ['firebrick', 'orangered', 'darkorange','gold', 'yellowgreen', 'lightgreen', 'turquoise', 'lightskyblue', 'plum', 'lightpink']

        plt.plot(iterations, scores, "-", linewidth=2, c=colors[n-1], label=n)
        # plt.show()

    def optimize_graph(self, n, iterations):
        """
        
        """
        plt.xlim(left=0, right=iterations)
        plt.ylim(top=0)
        plt.title("Hill-climber", fontweight='bold')
        plt.xlabel("Iterations", loc='right')
        plt.ylabel("Score", loc='top')
        plt.legend(range(1, n+1), title="Minimum amount of bonds changed per iteration", ncol=n//2)
        plt.show()


    def experiment(self, protein, iterations, sample_size=1, max_n=10):
        """
        
        """
        if self.prints: print("\n\nStarting score =", protein.score, "Length protein =", len(self.protein.sequence_list))
        
        start = time.time()
        if max_n > 10:
            print("ERROR: please only insert a max_n of 10 or smaller") #Make this a RaiseValueError
            return
        
        first_random_fold = copy.deepcopy(protein)
        best_protein = first_random_fold
        best_score = first_random_fold.score
        
        for n in range(1, max_n+1):
            if self.prints: print(f"\nStarting new N: {n}")

            all_scores = {}
            for iteration in range(iterations+1):
                all_scores[iteration] = []

            for sample in range(sample_size):
                protein_for_n, lowest_score_for_n, scores, improvement = self.run_n_iterations(first_random_fold, iterations, n)

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
        self.optimize_graph(max_n, iterations)
        if self.prints: print(f"\n\nBest score has become", best_score)
        if self.prints: print(f"Runtime experiment: {end-start} seconds.\n")
        best_protein.create_output()
        return best_protein
    