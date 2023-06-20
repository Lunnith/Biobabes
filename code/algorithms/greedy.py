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

    def __init__(self, protein, dimensions, splits = 3):
        self.protein = protein
        self.splits = splits
        self.dimensions = dimensions


        self.number_of_splits = (len(self.protein.sequence) - 1) // self.splits
        self.amino_left = (len(self.protein.sequence) - 1) - (self.number_of_splits * self.splits)
        self.used_coordinates_G = set()
        

        if self.dimensions == 2:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
        elif self.dimensions == 3:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))
    
    def all_bonds(self):

        # initiate i and keep in the while loop until i has gone trough all protein_sequence elements
        i = 0
        while i < len(self.protein.sequence):
        
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
        
        # initiate a set of coordinates that have been used within this direction state
        temp_list_coordinates = set()

        # set protein score to zero to determine the score of adding this part of the sequence
        self.protein.score = 0

        for index, direction in enumerate(state_directions):
        
            # for every element in the protein sequence, add and select the aminoacid
            self.protein.add_aminoacid(self.protein.sequence[i + index])
            self.acid = self.protein.sequence_list[-1]

            # create a bond using the acid, the previous acid and the direction
            self.protein.create_bond(self.acid, self.protein.sequence_list[-2], direction)

     
            # check if the acid location is not already occupied or a location that will be stuck
            if tuple(self.acid.location) not in self.used_coordinates_G and self.is_stuck(self.acid) == False and tuple(self.acid.location) not in temp_list_coordinates:
        
                # check the new interations and assign the acquired score to score direction
                self.acid.check_interactions(self.protein)
               
                # add the coordinates to the temporary list of coordinates
                temp_list_coordinates.add((tuple(self.acid.location)))
                
            else:

                # remove the amount of aminoacids that have been added and return False
                del self.protein.sequence_list[-(index + 1):]
                return False
 
        # remove the amount of aminoacids that have been added
        del self.protein.sequence_list[-(index + 1):]
        
        # return the protein score
        return self.protein.score

    def create_directions(self, size):
        
        directions = []
        for i in range(size):
            if self.dimensions == 2:
                directions.append([tuple((1, 0, 0, 1)), tuple((-1, 0, 0, -1)), tuple((0, 1, 0, 2)), tuple((0, -1, 0, -2))])
            else:
                directions.append([tuple((1, 0, 0, 1)), tuple((-1, 0, 0, -1)), tuple((0, 1, 0, 2)), tuple((0, -1, 0, -2)), tuple((0, 0, 1, 3)), tuple((0, 0, -1, -3))])
        self.list_directions = list(itertools.product(*directions))
        print(len(self.list_directions))

        return self.list_directions
        
    
    def all_bonds_kk(self):

        # set the first aminoacid at coordinate 0,0,0 and add to used coordinates
        self.protein.add_aminoacid(self.protein.sequence[0])
        self.acid = self.protein.sequence_list[0]
        self.acid.location = [0, 0, 0]
        self.used_coordinates_G.add((tuple(self.acid.location)))

        # initiate all directions
        self.list_directions = self.create_directions(self.splits)
        i = 1
        # iterate the lenght of the protein sequence (skipping the first) in steps of the split
        while i in range(len(self.protein.sequence)):
            print("ROUND", i)
            
            # initiate a valid random direction for the whole split in case non of the directions lead to a better score than 0
            best_state_directions = []

            # make a set for the coordinates used within this split to make sure this part of the protein does not go over itself
            self.used_coordinates_random = set()

            if i == 1 + (self.number_of_splits * self.splits):
                splits = self.amino_left
            else:
                splits = self.splits
            print('splits', splits)
            j = 0 
            while j in range(splits):
                print('j', j)
                # for every aminoacid, add to the protein and 
                self.protein.add_aminoacid(self.protein.sequence[i + j])
                self.acid = self.protein.sequence_list[-1]

                # append valid direction to the best state directions and add to the used random coordinates
                random_bond = self.random_bond()
                print(random_bond)
                if random_bond == False and j != 0:
                    self.used_coordinates_G.add(tuple(self.protein.sequence_list[-2].location))
                    del self.protein.sequence_list[-2:]
                    j -= 1
                    continue

                elif random_bond == False and j == 0:
                    i -= self.splits
                    self.used_coordinates_G.add(tuple(self.protein.sequence_list[-2].location))
                    del self.protein.sequence_list[-1]
                    break

                else:
                    best_state_directions.append(random_bond)
                    self.used_coordinates_random.add((tuple(self.acid.location)))
                    j += 1

            # delete the used aminoacids from the protein sequence list
            del self.protein.sequence_list[-(splits):]

            if random_bond == False:
                continue

            best_score = 0

            if splits == self.amino_left:
                list_directions = [self.directions]
            
            else:
                list_directions = self.list_directions

            for state_directions in list_directions:
                
                score = self.parts(state_directions, i)

                if score < best_score and score != False:
    
                    best_score = score
                    best_state_directions = state_directions
            print(best_score)
            for index, direction in enumerate(best_state_directions):

                self.protein.add_aminoacid(self.protein.sequence[i + index])
                self.acid = self.protein.sequence_list[-1]
                self.add_best_direction(direction)
            i += self.splits

    def add_best_direction(self, best_direction):

        # create a bond for the best direction for this acid and add the location to used coordinates
        self.protein.create_bond(self.acid, self.protein.sequence_list[-2], best_direction)
        #self.best_directions.append(best_direction)
        print(self.acid.location)
        self.used_coordinates_G.add(tuple(self.acid.location))

    def check_all_interactions(self):
        """
        Check all the interactions and adjust the score when an interaciton is found
        """
        # set the lists of bonds to empty again (because they were filled during directions)
        self.protein.hh_ch_bonds = []
        self.protein.cc_bonds = []
        self.protein.score = 0

        # go through the protein sequence list and check the interaction for every aminoacid
        for index, acid in enumerate(self.protein.sequence_list):
            acid.check_interactions(self.protein, index)
        
        # return the protein
        return self.protein


    def is_stuck(self, acid):
        """
        Check if an acid placement is stuck, does not have any direction to go in
        """
       
        # go trough all possible directions and add the direction to the acid location and check if it is in used coordinates
        for direction in self.directions:
            if tuple(list(map(add, acid.location, direction[0:3]))) not in self.used_coordinates_G:

                # if any of the directions is not occupied
                return False
    
        # return True because is_stuck is true
        return True
              

    def random_bond(self):
        """
        This function creates a random valid, random direction. If there is no valid direction, it will return False
        """
        #set best direction to a random direction and create the bond with this direction
        best_direction = random.choice(tuple(self.directions))
        self.protein.create_bond(self.acid, self.protein.sequence_list[-2], best_direction)

        # Until the acid location that is created is a valid location (not occupied or stuck) keep trying different directions
        tried_directions = set()
        while tuple(self.acid.location) in self.used_coordinates_G or self.is_stuck(self.acid) == True or tuple(self.acid.location) in self.used_coordinates_random:

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
        print(len(self.protein.sequence_list))



