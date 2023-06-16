# from code.classes.protein import Protein
from code.visualization.visualize import *
import random
import copy
import matplotlib.pyplot as plt
from operator import add

class Greedy():

    def __init__(self, protein, dimensions):
        self.protein = protein
        #self.splits = splits

        #self.best_directions = []
        self.used_coordinates_G = set()
        if dimensions == 2:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
        if dimensions == 3:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))
    
    def all_bonds(self):

        # initiate i and keep in the while loop until i has gone trough all protein_sequence elements
        i = 0
        while i < len(self.protein.sequence):
            print(i)
            # for every element in the protein sequence, add and select the aminoacid
            self.protein.add_aminoacid(self.protein.sequence[i])
            self.acid = self.protein.sequence_list[i]
            
            # add location for the first acid and add these to the used coordinates list
            if i == 0:
                self.acid.location = [0, 0, 0]
                self.used_coordinates_G.add((tuple(self.acid.location)))

            # for other aminoacids then the first create bond with previous acid
            elif i != 0:
                best_score = 0

                # create a random direction with random_bond to use if all the directions have the same score
                best_direction = self.random_bond()

                # check if there is no direction the acid can move in
                if best_direction == False:

                    # delete the acid of the current and previous iteration within the protein list
                    del self.protein.sequence_list[-2:]

                    ###del self.best_directions[-1]

                    # go back one iteration to not get stuck again
                    i -= 1

                    # skip the rest of the loop 
                    continue
                
                
                # go trough all directions to see which direction leads to the lowest score
                for direction in self.directions:
                    
                    # make a deepcopy of the protein and use this to index the acid
                    self.protein_temp = copy.deepcopy(self.protein)
                    self.acid_temp = self.protein_temp.sequence_list[i]
                    
                    # reset the score of the temporary protein to zero
                    self.protein_temp.score = 0
                    self.protein_temp.create_bond(self.acid_temp, self.protein_temp.sequence_list[i - 1], direction)

                    # check if the acid location is not already occupied or a location that will be stuck
                    if tuple(self.acid_temp.location) not in self.used_coordinates_G and self.is_stuck(self.acid_temp) == False:
                        
                        # check the new interations and assign the acquired score to score direction
                        self.acid_temp.check_interactions(self.protein_temp)
                        score_direction = self.protein_temp.score

                        # check if the score of the directon is smaller than the current best score
                        if score_direction < best_score:
                            
                            # change best score to the score of the direction and replace the best direction with the current direction
                            best_score = score_direction
                            best_direction = direction
                
                # add the best direction for the acid to the protein
                self.add_best_direction(best_direction)

            # go forward one iteration in the while loop
            i += 1

    def add_best_direction(self, best_direction):
        
        # create a bond for the best direction for this acid and add the location to used coordinates
        self.protein.create_bond(self.acid, self.protein.sequence_list[-2], best_direction)
        #self.best_directions.append(best_direction)
        self.used_coordinates_G.add(tuple(self.acid.location))

    def check_all_interactions(self):
        """
        Check all the interactions and adjust the score when an interaciton is found
        """
        for index, acid in enumerate(self.protein.sequence_list):
            acid.check_interactions(self.protein, index)
                
        return self.protein


    def is_stuck(self, acid):
        """
        Check if an acid is stuck, does not have any direction to go in
        """
        # counter initialization
        k = 0

        # go trough all possible directions and add the direction to the acid location and check if it is in used coordinates
        for direction in self.directions:
            if tuple(list(map(add, acid.location, direction[0:3]))) in self.used_coordinates_G:

                # if this direction is already occupied, add 1 to the counter
                k += 1

        # if all directions are occupied, add the acid.location to used coordinates because this location is stuck and therfore not usable
        if k == 4:
            self.used_coordinates_G.add(tuple(acid.location))

            # return True because is_stuck is true
            return True

        # if not all directions are occupied, return False
        return False
        

    def random_bond(self):
        """
        This function creates a random valid, random direction. If there is no valid direction, it will return False
        """
        #set best direction to a random direction and create the bond with this direction
        best_direction = random.choice(tuple(self.directions))
        self.protein.create_bond(self.acid, self.protein.sequence_list[-2], best_direction)

        # Until the acid location that is created is a valid location (not occupied or stuck) keep trying different directions
        tried_directions = set()
        while tuple(self.acid.location) in self.used_coordinates_G or self.is_stuck(self.acid) == True:

            # repeat trying new directions to find a valid one
            best_direction = random.choice(tuple(self.directions))
            self.protein.create_bond(self.acid, self.protein.sequence_list[-2], best_direction)
            tried_directions.add(best_direction)

            # if there all directions are checked and invalid, return False
            if tried_directions == self.directions:
                print('ERROR STUCK')
                return False
            
        # return a valid random direction 
        return best_direction
    
    def run(self):

        self.all_bonds()
        self.check_all_interactions()
        self.protein.create_output()
        
                



