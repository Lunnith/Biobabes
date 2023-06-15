from .randomise import random_assignment
import random

class Hill_climber():
    def __init__(self, protein, dimensions=3):
        self.protein = protein
        
        #Set options for dimensions
        if dimensions == 2:
            self.directions = {1: [1, 0, 0, 1], -1: [-1, 0, 0, -1], 2: [0, 1, 0, 2], -2: [0, -1, 0, -1]}
        if dimensions == 3:
            self.directions = {1: [1, 0, 0, 1], -1: [-1, 0, 0, -1], 2: [0, 1, 0, 2], -2: [0, -1, 0, -2], 3: [0, 0, 1, 3], -3: [0, 0, -1, -3]}


        #Initiate first random folding
        self.protein = random_assignment(protein, dimensions)
    
    def change_one_bond(self, protein, given_index=None):
        """
        To do: Remove old coordinate from protein.used_coordinates (which is a set) and add new coords.
                (did I accidentally fix this problem in the refold function?)
        To do: fix the testing if all directions have been tried (acid.location_valid check)
        To do: fix bug with the dictionary key = 0
        """
        if given_index == None:            
            index_changing_bond = random.randint(1, len(protein.sequence_list)-1)
        else:
            index_changing_bond = given_index
        
        print(index_changing_bond)

        location_previous_aminoacid = protein.sequence_list[index_changing_bond - 1].location
        old_location = protein.sequence_list[index_changing_bond].location
        old_step = protein.sequence_list[index_changing_bond].step
        print("Old step is", old_step)
        new_step = old_step
        acid = protein.sequence_list[index_changing_bond]

        #Get the old direction
        old_direction = self.directions[old_step]
        print("Old direction is", old_direction)
        

        #Determine new step
        tried_directions = set()
        # tried_directions.add(old_direction)
        acid.location_valid = False
        
        while acid.location_valid == False:
            while new_step == old_step:
                new_step = random.sample(self.directions.keys(), 1)[0]
            print("New step is", new_step)
            #Overwrite step of previous acid
            protein.sequence_list[index_changing_bond - 1].step = new_step
            #fix format of directions
            print("New directions of this aminoacid", self.directions[new_step])
            direction = tuple(self.directions[new_step])
            print("New direction as tuple:", direction)

            #Create new bond based on new direction
            protein.create_bond(acid, protein.sequence_list[index_changing_bond - 1], direction)
            new_location = protein.sequence_list[index_changing_bond].location
            # tried_directions.add(direction)

            # if protein can't fold anymore, return shorter folded protein
            #This doesnt work for dictionaries like this
            if tried_directions == self.directions:
                print('ended')
                return protein, False     

            print(location_previous_aminoacid, old_location, index_changing_bond, new_location)
            return protein, index_changing_bond
    
    def refold(self, protein, index):
        """
        Keep all directions of the protein and fold it exactly the same way, but with the one changed bond.
        """
        if index != False:
            #Reset hh/ch/cc bonds
            #Maybe put this in the change_one_bond already?
            protein.hh_ch_bonds = []
            protein.cc_bonds = []

            used_coords = set()
            #Retrieve folded protein up until the changed bond
            for acid in range(0, index):
                print(tuple(protein.sequence_list[acid].location))
                used_coords.add(tuple(protein.sequence_list[acid].location))
            print(used_coords)
            #Refold from the changed bond
            for acid in range(index-1, len(protein.sequence_list)):
                direction = self.directions[protein.sequence_list[acid-1].step]
                protein.create_bond(protein.sequence_list[acid], protein.sequence_list[acid-1], direction)
                used_coords.add(tuple(protein.sequence_list[acid].location))
            print(used_coords)
            protein.used_coordinates = used_coords


            
            #check_validity()
            return protein
        else:
            return False

    def check_score(self, protein):
        """
        Calculate the score of this new fold
        """
        #for acid in protein.sequence_list:
            #acid.check_interactions()
        pass

    def check_validity(self, protein):
        """
        Check if the protein has folded over itself
        """
        # return True/False, double_coords (list) (if True, list is empty)
        pass

    def change_n_bonds(self, protein, n):
        """
        Loop over all functions n times
        while check_validity gives >0, first change these bonds
        """
        for i in range(n):
            if len(double_coords) == 0:
                protein, index_changed_bond = self.change_one_bond(protein)
                self.refold(protein, index_changed_bond)
                #self.check_validity
            else:
                #Change the double coords in order from first to last
                self.change_one_bond(protein, given_index=double_coords[0])
        #If after n changes, there are still double bonds, this initial change in bond is invalid
        #Else, check_score()
        #return protein, score



