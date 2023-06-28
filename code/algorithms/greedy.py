from code.classes.protein import Protein
from code.classes.aminoacid import Aminoacid
import random
from operator import add
import itertools

class Greedy():
    """
    A class to do a greedy algorithm and beam search depending on the input. 

    Attributes:
    -----------
    Protein: class object
        The protein with a sequence to initiate the acids and to determine/make the best folding protein
    
    Acid: class object
        Aminoacid that is added to the protein
    
    Splits: integer
        Splits determine the width of the beam in a beam search. When set to 1, it will act as a normal greedy algorithm
    
    Dimensions: integer
        Dimensions can be set to 2 or 3 for 2D and 3D folding respectively. Dimensions determine the initiation of the set of directions
    
    Directions: set
        Dependend on the dimensions, the set of directions contains 4 or 6 possible directions with the change in x, y or z and the step for the output
    
    Amino before: integer
        Determines after how many aminoacids the beam search starts, this is to shift the beam search so it is not the same every time
    
    Number of splits: integer
        This is the number of times a beam fits in the sequence dependend on the amount of aminoacids before and the sequence lenght
    
    Amino left: integer
        This is the amount of aminoacids left after the beam search, determined by the number of splits and aminoacids before

    Used coordinates G: set
        This is a list of coordinates that are occupied or lead to a location that is stuck

    list directions: list of lists
        There are 3 directions lists with all possible states of a sequence part, one for the amount of aminoacids of before beam search, 
        one for the splits and one after

    Methods:
    -----------
    Create directions(size):
        
    All bonds(): 
        Determine the best directions for a sequence part leading to lowest score

    Parts(state_directions, i):
        Determine the score of adding directions for one state

    Add best direction(best_direction, place):
        Adds the best direction for a sequence part to the protein

    check all interactions():
        Checks interactions and add the scores of every aminoacid with the previous acids

    is stuck(acid):
        Checks if an acid is stuck (surrounded by occupied coordinates)

    run:
        Runs the algorithm after initiation
    """
    def __init__(self, protein: Protein, dimensions: int, splits: int = 1, before: int = 0) -> None:

        # define protein, splits and dimensions
        self.protein = protein
        self.splits = splits
        self.dimensions = dimensions

        self.states = 0

        # set the amount of aminoacids before the sequence parts start (0 means it start immediatly)
        self.amino_before = before

        # calculate how much aminoacids er are left after deviding the sequence parts
        self.number_of_splits = (len(self.protein.sequence) - 1 - self.amino_before) // self.splits
        self.amino_left = (len(self.protein.sequence) - 1) - (self.number_of_splits * self.splits) - self.amino_before

        # initiate een set voor de gebruikte coordinaten
        self.used_coordinates_G = set()

        # based on the dimensions, give self directions 4 (x and y) or 6 (x, y and z) possible directions
        if self.dimensions == 2:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
        elif self.dimensions == 3:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))

        # initiate all directions for the splits,the aminoacids before and for the leftover aminoacids at the end
        self.list_directions = self.create_directions(self.splits)
        self.list_directions_left = self.create_directions(self.amino_left)
        self.list_directions_before = self.create_directions(self.amino_before)


    def create_directions(self, size: int) -> list[tuple]:
        """
        This function creates all direction combinations for every state of a sequence part
        """
        # create an empty directions list and fill it with all possible directions, times the amount of the size of the sequence part
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
        
    
    def all_bonds(self) -> None:
        """
        This function goes trough all the states and selects one of the states with the lowest score
        """
        # set the first aminoacid at coordinate 0,0,0 and add to used coordinates
        self.protein.add_aminoacid(self.protein.sequence[0])
        self.acid = self.protein.sequence_list[0]
        self.acid.location = [0, 0, 0]
        self.used_coordinates_G.add((tuple(self.acid.location)))

        # iterate the lenght of the protein sequence (skipping the first) in steps of the split
        i = 1
        while i in range(len(self.protein.sequence)):
            
            # initiate an empty list for best state directions with the same (best) score ## weg?
            best_state_directions = []

            # initiate an empty list for states with the same (best) score
            best_states = []

            # make a set for the coordinates used within this split to make sure this part of the protein does not go over itself
            self.used_coordinates_random = set()

            # Depending on where in the sequence we are, set splits and directions list to the accompanying the right amount of aminoacids
            if i == 1 + (self.number_of_splits * self.splits) + self.amino_before:
                splits = self.amino_left
                list_directions = self.list_directions_left
            elif i == 1 and self.amino_before > 0:
                splits = self.amino_before
                list_directions = self.list_directions_before
            else:
                splits = self.splits
                list_directions = self.list_directions

            # go through all states and compute the score
            best_score = 1
            for state_directions in list_directions:

                # add one state to the states counter
                self.states += 1
                # compute the score using the parts function
                score = self.parts(state_directions, i)
                if score == best_score and score != None:
                    best_states.append(state_directions)        
                
                # determine the best score and best direction
                elif score != None and score < best_score:
                    best_states = []
                    best_states.append(state_directions)
                    best_score = score

            # when the best score did not change, this means every score was false and the protein is stuck
            if best_score == 1:

                # go back one split and skip the rest of this iteration
                i -= self.splits
                del self.protein.sequence_list[-(self.splits):]
                continue
            
            # take one of the states with the lowest score
            best_state_directions = random.choice(best_states)

            # add the aminoacids to self.protein with the directions of the best score
            for index, direction in enumerate(best_state_directions):
                self.add_best_direction(direction, place = (i + index))
            
            # go to the next sequence part
            i += splits



    def parts(self, state_directions: tuple(tuple()), i: int) -> int:
        """
        This function makes the bonds for a state of a sequence part and returns the score
        """
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
    


    def add_best_direction(self, best_direction: tuple, place: int) -> None:
        """"
        This function creates a bond for the best direction for an acid and adds the location to used coordinates 
        """
        self.protein.add_aminoacid(self.protein.sequence[place])
        self.acid = self.protein.sequence_list[-1]
        self.protein.create_bond(self.acid, self.protein.sequence_list[-2], best_direction)
        self.used_coordinates_G.add(tuple(self.acid.location))



    def check_all_interactions(self) -> Protein:
        """
        Check all the interactions and adjust the score when an interaciton is found
        """
        # set the lists of bonds to empty again and score to zero (because they were filled during directions)
        self.protein.hh_ch_bonds = []
        self.protein.cc_bonds = []
        self.protein.score = 0

        # go through the protein sequence list and check the interaction for every aminoacid
        for index, acid in enumerate(self.protein.sequence_list):
            acid.check_interactions(self.protein, index)
        
        # return the protein
        return self.protein


    def is_stuck(self, acid: Aminoacid) -> bool:
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
    
    def run(self) -> None:
        """
        This function runs the script by first setting all the (best) bonds, then checking the interactions and creating the output
        """
        self.all_bonds()
        self.check_all_interactions()
        #self.protein.create_output()



