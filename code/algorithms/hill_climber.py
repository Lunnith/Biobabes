from .randomise import random_assignment
import random
import time

class Hill_climber():
    def __init__(self, protein, dimensions=3):
        start = time.time()
        self.protein = protein
        self.double_coords = {}
        
        #Set options for dimensions
        if dimensions == 2:
            self.directions = {1: [1, 0, 0, 1], -1: [-1, 0, 0, -1], 2: [0, 1, 0, 2], -2: [0, -1, 0, -1]}
        if dimensions == 3:
            self.directions = {1: [1, 0, 0, 1], -1: [-1, 0, 0, -1], 2: [0, 1, 0, 2], -2: [0, -1, 0, -2], 3: [0, 0, 1, 3], -3: [0, 0, -1, -3]}


        #Initiate first random folding
        self.protein = random_assignment(protein, dimensions)
        self.lowest_score = protein.score
        self.best_fold = protein
        end = time.time()
        print(f"Runtime __init__: {end-start} seconds.")
    
    def change_one_bond(self, protein, given_index=None):
        """
        To do: Remove old coordinate from protein.used_coordinates (which is a set) and add new coords.
                (did I accidentally fix this problem in the refold function?)
        To do: fix the testing if all directions have been tried (acid.location_valid check)
        To do: fix bug with the dictionary key = 0
        """
        start = time.time()
        if given_index == None:            
            index_changing_bond = random.randint(1, len(protein.sequence_list)-1)
        else:
            index_changing_bond = given_index

        location_previous_aminoacid = protein.sequence_list[index_changing_bond - 1].location
        old_location = protein.sequence_list[index_changing_bond].location
        old_step = protein.sequence_list[index_changing_bond].step
        new_step = old_step
        acid = protein.sequence_list[index_changing_bond]

        #Get the old direction
        old_direction = self.directions[old_step]

        #Determine new step
        tried_directions = set()
        # tried_directions.add(old_direction)
        acid.location_valid = False
        
        while acid.location_valid == False:
            while new_step == old_step:
                new_step = random.sample(self.directions.keys(), 1)[0]
            #Overwrite step of previous acid
            protein.sequence_list[index_changing_bond - 1].step = new_step
            #fix format of directions
            direction = tuple(self.directions[new_step])

            #Create new bond based on new direction
            protein.create_bond(acid, protein.sequence_list[index_changing_bond - 1], direction)
            new_location = protein.sequence_list[index_changing_bond].location
            # tried_directions.add(direction)

            # if protein can't fold anymore, return shorter folded protein
            #This doesnt work for dictionaries like this
            if tried_directions == self.directions:
                print('ended')
                return protein, False     

            end = time.time()
            print(f"Runtime change_one_bond: {end-start} seconds.")
            return protein, index_changing_bond
    
    def refold(self, protein, index):
        """
        Keep all directions of the protein and fold it exactly the same way, but with the one changed bond.
        """
        start = time.time()
        if index != False:
            #Reset hh/ch/cc bonds
            #Maybe put this in the change_one_bond already (or in change_n_bonds)?
            protein.hh_ch_bonds = []
            protein.cc_bonds = []

            used_coords = set()
            #Retrieve folded protein up until the changed bond
            for acid in range(0, index):
                used_coords.add(tuple(protein.sequence_list[acid].location))
            #Refold from the changed bond
            for acid in range(index-1, len(protein.sequence_list)):
                direction = self.directions[protein.sequence_list[acid-1].step]
                protein.create_bond(protein.sequence_list[acid], protein.sequence_list[acid-1], direction)
                used_coords.add(tuple(protein.sequence_list[acid].location))
            protein.used_coordinates = used_coords


            
            #check_validity()
            end = time.time()
            print(f"Runtime Refold: {end-start} seconds.")
            return protein

    def check_score(self, protein):
        """
        Calculate the score of this new fold
        YET TO FIX BUG: Also counts the normal bonds in the score
            Problem found: Check_interactions is built thinking that future acids are not yet present.
            To do: find a way to give check_interactions only the protein up untill the acid of interest.
        """
        protein.score = 0
        start = time.time()
        for acid in protein.sequence_list:
            acid.check_interactions(protein)

        end = time.time()
        print(f"Runtime check_score: {end-start} seconds.")

        

    def check_validity(self, protein):
        """
        Check if the protein has folded over itself
        """
        start = time.time()
        self.double_coords = {}
        used_coords = set()
        for acid in protein.sequence_list:
            if tuple(acid.location) not in used_coords:
                used_coords.add(tuple(acid.location))
            else:
                index = protein.sequence_list.index(acid)
                self.double_coords[index] = tuple(acid.location)
        # return self.double_coords (dict, keys = index, values = coords)
        end = time.time()
        print(f"Runtime check_validity: {end-start} seconds")

        return self.double_coords

    def change_n_bonds(self, protein, n):
        """
        Loop over all functions n times
        while check_validity gives >0, first change these bonds

        #NOTE TO SELF: zorg dat ie alleen de eerste dubbele vouwing fixt, pas als die goed is de rest
            nu gebeurt het vaak dat ie binding 41, 48 en 51 fixt, 
            om vervolgens te zien dat 41 nogsteeds niet correct is, dus die opnieuw aan te passen,
            Waardoor uiteindelijk 51 weer dubbel zit.
        """
        start = time.time()
        for i in range(n):
            protein, index_changed_bond = self.change_one_bond(protein)
            self.refold(protein, index_changed_bond)

        #Fix the aminoacids that have been folded incorrectly
        while len(self.check_validity(protein)) > 0:
            print("Indexes of double acids =", self.double_coords.keys())
            for acid in self.double_coords.keys():
                protein, index_changed_bond = self.change_one_bond(protein, given_index=acid)
                self.refold(protein, index_changed_bond)
            # TO DO: Create safety to not get stuck.

        end = time.time()
        print(f"Runtime change_n_bonds: {end-start} seconds.")
        return protein



    def run_n_iterations(self, protein, iterations, bonds):
        """
        runs the change_n_bonds
        """
        print("Starting score =", self.lowest_score)
        print("\n")
        for n in range(iterations):
            new_protein = self.change_n_bonds(protein, bonds)
            self.check_score(new_protein)

            print("\n Score after n changed bonds", protein.score)
            if protein.score < self.lowest_score:
                print("Score updated to", protein.score, "\n\n")
                self.lowest_score = protein.score
                self.best_fold = protein

        return self.best_fold, self.lowest_score