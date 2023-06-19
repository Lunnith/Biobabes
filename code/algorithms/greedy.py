# from code.classes.protein import Protein
from code.visualization.visualize import *
from code.classes.protein import Protein
import random
import copy
from code.algorithms.depth_first import DepthFirst
import matplotlib.pyplot as plt
from operator import add
import itertools

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

    def parts(self, state_directions, i):
        
        remove_k = 0
        sequence_placement = i
        temp_list_coordinates = []

        for direction in state_directions:
            remove_k += 1
            #print("sequence placement", sequence_placement)
            # for every element in the protein sequence, add and select the aminoacid
            #print("KJKJKJKJ", self.protein.sequence[sequence_placement])
            self.protein.add_aminoacid(self.protein.sequence[sequence_placement])
            self.acid = self.protein.sequence_list[-1]
            
            # make a deepcopy of the protein and use this to index the acid
            self.protein_temp = copy.deepcopy(self.protein)
            self.acid_temp = self.protein_temp.sequence_list[-1]
            
            # reset the score of the temporary protein to zero
            
            #print("kkkk",self.protein_temp.sequence_list[-2].location)
           # print(direction)
            self.protein_temp.create_bond(self.acid_temp, self.protein_temp.sequence_list[-2], direction)
            self.protein.create_bond(self.acid, self.protein.sequence_list[-2], direction)
            #print(self.acid_temp.location)

            # check if the acid location is not already occupied or a location that will be stuck
            if tuple(self.acid_temp.location) not in self.used_coordinates_G and self.is_stuck(self.acid_temp) == False and self.acid_temp.location not in temp_list_coordinates:
                
                # check the new interations and assign the acquired score to score direction
                self.acid_temp.check_interactions(self.protein_temp)
                temp_list_coordinates.append(self.acid_temp.location)
                
            else:
                del self.protein.sequence_list[-remove_k:]
                return False
            
            sequence_placement += 1

        del self.protein_temp.sequence_list[-remove_k:]
        del self.protein.sequence_list[-remove_k:]
        print(len(self.protein.sequence_list))
        score = self.protein_temp.score
        self.protein_temp.score = 0
        return score

    def create_directions(self):
        
        directions = []
        for i in range(3):
            directions.append([tuple((1, 0, 0, 1)), tuple((-1, 0, 0, -1)), tuple((0, 1, 0, 2)), tuple((0, -1, 0, -2))])
        self.list_directions = list(itertools.product(*directions))
    
    def all_bonds_kk(self):
        self.protein.add_aminoacid(self.protein.sequence[0])
        self.acid = self.protein.sequence_list[0]
        self.acid.location = [0, 0, 0]
        self.used_coordinates_G.add((tuple(self.acid.location)))
        self.create_directions()
        
        for i in range(1, len(self.protein.sequence), 3):
            print("ROUND", i)
            #part_sequence = self.protein.sequence[i:i+3]
            best_state_directions = random.choice(self.list_directions) ###MOET VALID OPLOSSING ZIJN
            #self.acid_steps = {1: tuple((1,0,0,1)), -1: tuple((-1, 0, 0, -1)), 2: tuple((0, 1, 0, 2)) , -2: tuple((0, -1, 0, -2))}
            for state_directions in self.list_directions:
                best_score = 0
                score = self.parts(state_directions, i)
                #print("SCORE", score)
                if score == False:
                    continue
                
                else:
                    if score < best_score:
                    
                        best_score = score
                        best_state_directions = state_directions
            counter_score = i
            for direction in best_state_directions:
                #print("best_state_directions", best_state_directions)
                self.protein.add_aminoacid(self.protein.sequence[counter_score])
                self.acid = self.protein.sequence_list[-1]
                self.add_best_direction(direction)
                counter_score += 1
        

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
    
    def run_k(self):
        self.all_bonds_kk()
        self.check_all_interactions()
        self.protein.create_output()



