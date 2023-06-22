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

    def __init__(self, protein, dimensions, splits = 3, before = 0):

        # define protein, splits and dimensions
        self.protein = protein
        self.splits = splits
        self.dimensions = dimensions

        # calculate how much aminoacids er over blijven als je het in stukjes hebt verdeeld
        self.amino_before = before
        self.number_of_splits = (len(self.protein.sequence) - 1 - self.amino_before) // self.splits
        self.amino_left = (len(self.protein.sequence) - 1) - (self.number_of_splits * self.splits) - self.amino_before

        # initiate een set voor de gebruikte coordinaten
        self.used_coordinates_G = set()

        # based on the dimensions, give self directions 4 (x and y) or 6 (x, y and z) possible directions
        if self.dimensions == 2:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
        elif self.dimensions == 3:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))

        # initiate all directions for the parts and for the left over aminoacids
        self.list_directions = self.create_directions(self.splits)
        self.list_directions_left = self.create_directions(self.amino_left)
        self.list_directions_before = self.create_directions(self.amino_before)


    def create_directions(self, size):
        """
        This function creates all direction combinations for every state of a sequence part
        """
        #create an empty directions list and fill it with all possible directions times the amount of the size of the sequence part
        directions = []
        for i in range(size):
            if self.dimensions == 2:
                directions.append([tuple((1, 0, 0, 1)), tuple((-1, 0, 0, -1)), tuple((0, 1, 0, 2)), tuple((0, -1, 0, -2))])
            else:
                directions.append([tuple((1, 0, 0, 1)), tuple((-1, 0, 0, -1)), tuple((0, 1, 0, 2)), tuple((0, -1, 0, -2)), tuple((0, 0, 1, 3)), tuple((0, 0, -1, -3))]) 
        
        # create all possible combinations of the lists in directions
        list_directions = list(itertools.product(*directions))

        # return the list with all directions
        return list_directions
        
    
    def all_bonds(self):

        # set the first aminoacid at coordinate 0,0,0 and add to used coordinates
        self.protein.add_aminoacid(self.protein.sequence[0])
        self.acid = self.protein.sequence_list[0]
        self.acid.location = [0, 0, 0]
        self.used_coordinates_G.add((tuple(self.acid.location)))

        # iterate the lenght of the protein sequence (skipping the first) in steps of the split
        i = 1
        while i in range(len(self.protein.sequence)):
            print("ROUND", i)
            best_states = []
            
            # initiate a valid random direction for the whole split in case non of the directions lead to a better score than 0
            best_state_directions = []

            # make a set for the coordinates used within this split to make sure this part of the protein does not go over itself
            self.used_coordinates_random = set()

            # if it is the round of left aminoacids, set splits to the amount of left aminoacids and list_directions to directions of amino left
            if i == 1 + (self.number_of_splits * self.splits) + self.amino_before:
                splits = self.amino_left
                list_directions = self.list_directions_left
            elif i == 1 and self.amino_before > 0:
                splits = self.amino_before
                list_directions = self.list_directions_before
            else:
                splits = self.splits
                list_directions = self.list_directions

            
            # for every aminoacid in a split, create a bond
            # j = 0 
            # while j in range(splits):

    
            #     # for every aminoacid, add to the protein and select the acid
            #     self.protein.add_aminoacid(self.protein.sequence[i + j])
            #     self.acid = self.protein.sequence_list[-1]

            #     # append valid direction to the best state directions and add to the used random coordinates
            #     #random_bond = self.random_bond(j)

            #     for direction in self.directions:
            #         if tuple(list(map(add, self.protein.sequence_list[-2].location, direction[0:3]))) not in self.used_coordinates_G:
            #             random_bond = True
            #             print(self.protein.sequence_list[-2].location)
            #         else:
            #             if j == 0:
            #                 random_bond = '0'
            #             else:
            #                 random_bond = '1'
            #             print(self.protein.sequence_list[-2].location)



            #     # if there is no bond possible within the first aminoacid of a sequence part, go back one sequence part and break this loop
            #     if random_bond == '0':
            #         i -= splits
            #         break
                
            #     # if there is no bond possible after the first aminoacid of a sequence, go back one iteration in the sequence
            #     elif random_bond == '1':
            #         j -= 1

            #     # add the random direction to the best state directions and add to the temporary used coordinates in random
            #     else:
            #         #best_state_directions.append(random_bond)
            #         #self.used_coordinates_random.add((tuple(self.acid.location)))

            #         # go to the next iteration in the loop
            #         j += 1

            # delete the added aminoacids from the protein sequence list
            #del self.protein.sequence_list[-(splits):]

            # if there is no possible direction, skip the rest of this iteration
            # if random_bond == '0' or random_bond == '1':
            #     continue

            # go through all states and compute the score
            best_score = 1
            for state_directions in list_directions:
                
                # compute the score using the parts function
                score = self.parts(state_directions, i)
                if score == best_score and score != None:
                    best_states.append(state_directions)        
                
                # determine the best score and best direction
                elif score != None and score < best_score:
                    best_states = []
                    best_states.append(state_directions)
                    best_score = score
            print(best_score)
            if best_score == 1:
                i -= self.splits
                del self.protein.sequence_list[-(self.splits):]
                continue

            best_state_directions = random.choice(best_states)

            # add the aminoacids to self.protein with the directions of the best score
            for index, direction in enumerate(best_state_directions):
                self.add_best_direction(direction, place = (i + index))
            
            
            # go to the next sequence part
            i += splits



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
        
                # check the new interations and add the location to temporary list of coordinates
                self.acid.check_interactions(self.protein)
                temp_list_coordinates.add((tuple(self.acid.location)))

            else:
                # remove the amount of aminoacids that have been added and return False
                del self.protein.sequence_list[-(index + 1):]
                return None
 
        # remove the amount of aminoacids that have been added
        del self.protein.sequence_list[-(index + 1):]
        
        # return the protein score
        return self.protein.score
    


    def add_best_direction(self, best_direction, place):

        # create a bond for the best direction for this acid and add the location to used coordinates
        self.protein.add_aminoacid(self.protein.sequence[place])
        self.acid = self.protein.sequence_list[-1]
        self.protein.create_bond(self.acid, self.protein.sequence_list[-2], best_direction)
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
              

    def random_bond(self, j):
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
                if j == 0:

                    # add the coordinate of the stuck location to the used coordinates and remove last added amioacid
                    self.used_coordinates_G.add(tuple(self.protein.sequence_list[-2].location))
                    del self.protein.sequence_list[-1]

                    # return '0'
                    return '0'
                
                # if there is no bond possible after the first aminoacid of a sequence, go back one iteration in the sequence
                elif j != 0:

                    # add the coordinate of the stuck location to the used coordinates and remove last and previous added amioacid
                    self.used_coordinates_G.add(tuple(self.protein.sequence_list[-2].location))
                    del self.protein.sequence_list[-2:]

                    # return '1'
                    return '1'
            
        # return a valid random direction 
        return best_direction
    
    def run(self):
        self.all_bonds()
        self.check_all_interactions()
        self.protein.create_output()
        



