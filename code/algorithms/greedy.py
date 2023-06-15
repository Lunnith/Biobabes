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
       # if dimensions == 2:
        self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
        # if dimensions == 3:
        #     self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))
    
    def all_bonds(self):

        for i in range(len(self.protein.sequence)):

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

                # go trough all directions to see which direction leads to the lowest score
                for direction in self.directions:
                  
                    self.protein_temp = copy.deepcopy(self.protein)
                    self.acid_temp = copy.deepcopy(self.acid)
                    self.protein_temp.score = 0
                    self.protein_temp.create_bond(self.acid_temp, self.protein_temp.sequence_list[i - 1], direction)
                
                    if tuple(self.acid_temp.location) in self.used_coordinates_G:
                        pass
                    else:
                        
                        self.acid_temp.check_interactions(self.protein_temp)
                        score_direction = self.protein_temp.score
                        print(i)
                        #print('score_direct', score_direction)

                        if score_direction < best_score:
                            print('inside', i)
                            best_score = score_direction
                            best_direction = direction
                            print('best_score', best_score)
                            print('best_direction', best_direction)

               
                self.protein.create_bond(self.acid, self.protein.sequence_list[i - 1], best_direction)       
                self.best_directions.append(best_direction)
                self.acid.check_interactions(self.protein)
                self.used_coordinates_G.add(tuple(self.acid.location))
                    #print(self.used_coordinates_G)
            
                               
    
    def random_bond(self, i):

        best_direction = random.choice(tuple(self.directions))
        self.protein.create_bond(self.acid, self.protein.sequence_list[i - 1], best_direction)

        tried_directions = set()
        while tuple(self.acid.location) in self.used_coordinates_G:

            best_direction = random.choice(tuple(self.directions))
            self.protein.create_bond(self.acid, self.protein.sequence_list[i - 1], best_direction)
            tried_directions.add(best_direction)

            # if protein can't fold anymore, return shorter folded protein
            if tried_directions == self.directions:
                return self.protein
        
        return best_direction
                




# def greedy(protein):
    
#     best_score = 0
#     best_directions = [] 
#    # if dimensions == 2:
#     directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
#     # if dimensions == 3:
#     #     directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))
#     used_coordinates_greedy = set()
#     for i in range(len(protein.sequence)):
#         protein.add_aminoacid(protein.sequence[i])
#         acid = protein.sequence_list[i]
        
#         # location of first aminoacid is (0,0,0)
#         if i == 0:
#             acid.location = [0,0,0]
        
#            # protein.used_coordinates.add((tuple(acid.location)))
#             used_coordinates_greedy.add((tuple(acid.location)))
           

#         # for other aminoacids then the first create bond with previous acid
#         elif i != 0:

#             best_score = 0
#             # set best direction to a random choice (still have to make it a valid one)
            
#             best_direction = random.choice(tuple(directions))
#             protein.create_bond(acid, protein.sequence_list[i - 1], best_direction)
#             tried_directions = set()
#             while tuple(acid.location) in used_coordinates_greedy:

#                 best_direction = random.choice(tuple(directions))
#                 protein.create_bond(acid, protein.sequence_list[i - 1], best_direction)
#                 tried_directions.add(best_direction)

#                 # if protein can't fold anymore, return shorter folded protein
#                 if tried_directions == directions:
#                     return protein
           

#             # go trough all directions to see which direction leads to the lowest score
#             for direction in directions:

#                 protein.create_bond(acid, protein.sequence_list[i - 1], direction)
#                 print(acid.location_valid)
#                 if tuple(acid.location) in used_coordinates_greedy:
#                 #if acid.location_valid == False:
#                     pass
#                 else:
#                     protein_temp = copy.deepcopy(protein)
#                     acid_temp = copy.deepcopy(acid)
#                     acid_temp.check_interactions(protein_temp)
#                     score_direction = protein_temp.score

#                     if score_direction < best_score:
#                         print("inside")
#                         best_score = score_direction
#                         best_direction = direction
                        
#                     print('location_dir', acid.location)
                
#             protein.create_bond(acid, protein.sequence_list[i - 1], best_direction)       
#             best_directions.append(best_direction)
#             acid.check_interactions(protein)
#             print('location', acid.location)
#             coordinates = tuple(acid.location)
#             used_coordinates_greedy.add((tuple(acid.location)))
           

#         print('used_cor', used_coordinates_greedy)
#     return best_directions


