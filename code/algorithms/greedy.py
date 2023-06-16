# from code.classes.protein import Protein
# from code.visualization.visualize import *
import random
import copy
import matplotlib.pyplot as plt
from operator import add

class Greedy():

    def __init__(self, protein):
        self.protein = protein
        #self.splits = splits

        self.best_directions = []
        self.used_coordinates_G = set()
        #dimensions = 2
        #if dimensions == 2:
        self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
        # if dimensions == 3:
        #     self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))
    
    def all_bonds(self):
        i = 0
        while i < len(self.protein.sequence):
        #for i in range(len(self.protein.sequence)):
            print(i)
            # for every element in the protein sequence, add an aminoacid
            self.protein.add_aminoacid(self.protein.sequence[i])
            self.acid = self.protein.sequence_list[i]
            

            if i == 0:
                self.acid.location = [0, 0, 0]
                self.used_coordinates_G.add((tuple(self.acid.location)))
            

            # for other aminoacids then the first create bond with previous acid
            elif i != 0:
                best_score = 0
                best_direction = self.random_bond(i)
                if best_direction == False:
                    del self.protein.sequence_list[-2:]
                    
                    i = i-1
                    
                    #self.used_coordinates_G.add(tuple(acid.location))
                    continue
                
                
                # go trough all directions to see which direction leads to the lowest score
                for direction in self.directions:
                    
                    self.protein_temp = copy.deepcopy(self.protein)
                    self.acid_temp = self.protein_temp.sequence_list[i]
                
                    self.protein_temp.score = 0
                    self.protein_temp.create_bond(self.acid_temp, self.protein_temp.sequence_list[i - 1], direction)
                
                    if tuple(self.acid_temp.location) in self.used_coordinates_G or self.check_if_stuck(self.acid_temp) == False:
                        
                        pass
                    else:
                
                        self.acid_temp.check_interactions(self.protein_temp)
                        score_direction = self.protein_temp.score

                        if score_direction < best_score:
                        
                            best_score = score_direction
                            best_direction = direction

                self.add_best_direction(best_direction, i)
            i = i + 1
                               
    def add_best_direction(self, best_direction, i):

        self.protein.create_bond(self.acid, self.protein.sequence_list[i - 1], best_direction)
        self.best_directions.append(best_direction)
        self.acid.check_interactions(self.protein)
        self.used_coordinates_G.add(tuple(self.acid.location))

    def check_if_stuck(self, amino):
        k = 0
        for direction in self.directions:
            if tuple(list(map(add, amino.location, direction[0:3]))) in self.used_coordinates_G:
                k += 1
        if k == 4:
            self.used_coordinates_G.add(tuple(amino.location))
            return False
        else:
            return True
        

    def random_bond(self, i):

        best_direction = random.choice(tuple(self.directions))
        self.protein.create_bond(self.acid, self.protein.sequence_list[i - 1], best_direction)

        tried_directions = set()
        while tuple(self.acid.location) in self.used_coordinates_G or self.check_if_stuck(self.acid) == False:

            best_direction = random.choice(tuple(self.directions))
            self.protein.create_bond(self.acid, self.protein.sequence_list[i - 1], best_direction)
            tried_directions.add(best_direction)

            # if protein can't fold anymore, return shorter folded protein
            if tried_directions == self.directions:
                print('ERROR STUCK')

                return False
        
        return best_direction
                



