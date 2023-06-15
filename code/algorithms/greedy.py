# from code.classes.protein import Protein
# from code.visualization.visualize import *
import random
import copy
import matplotlib.pyplot as plt
from operator import add

# class Greedy():

#     def __init__(self, protein, splits):
#         self.protein = protein
#         self.splits = splits
#         self.best_score = 0
#         self.best_directions = []
#         # if dimensions == 2:
#         directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
#         # if dimensions == 3:
#         #     directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))
#     def all_bonds(self):

#         for i in range(len(self.protein.sequence)):

#             self.protein.add_aminoacid(self.protein.sequence[i])

#             acid = protein.sequence_list[i]
#             if i == 0:
#                 acid.location = [0, 0, 0]
            
#                 protein.used_coordinates.add((tuple(acid.location)))
            

#             # for other aminoacids then the first create bond with previous acid
#             elif i != 0:

#                 # set best direction to a random choice (still have to make it a valid one)
#                 best_direction = random.choice(tuple(directions))

#                 # go trough all directions to see which direction leads to the lowest score
#                 for direction in directions:


def greedy(protein):
    
    best_score = 0
    best_directions = [] 
   # if dimensions == 2:
    directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
    # if dimensions == 3:
    #     directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))
    used_coordinates_greedy = set()
    for i in range(len(protein.sequence)):
        protein.add_aminoacid(protein.sequence[i])
        acid = protein.sequence_list[i]
        
        # location of first aminoacid is (0,0,0)
        if i == 0:
            acid.location = [0, 0, 0]
        
            protein.used_coordinates.add((tuple(acid.location)))
           

        # for other aminoacids then the first create bond with previous acid
        elif i != 0:
            print(i)
          #  tried_directions = set()
            best_score = 0
            # set best direction to a random choice (still have to make it a valid one)
            
            best_direction = random.choice(tuple(directions))
            protein.create_bond(acid, protein.sequence_list[i - 1], best_direction)
            tried_directions = set()
            while tuple(acid.location) in used_coordinates_greedy:
            #while acid.location_valid == False:

                best_direction = random.choice(tuple(directions))
                protein.create_bond(acid, protein.sequence_list[i - 1], best_direction)
                tried_directions.add(best_direction)

                # if protein can't fold anymore, return shorter folded protein
                if tried_directions == directions:
                    return protein
            #protein.used_coordinates.remove((tuple(acid.location)))

            # go trough all directions to see which direction leads to the lowest score
            for direction in directions:

                protein.create_bond(acid, protein.sequence_list[i - 1], direction)
                print(acid.location_valid)
                if tuple(acid.location) in used_coordinates_greedy:
                #if acid.location_valid == False:
                    pass
                else:
                    protein_temp = copy.deepcopy(protein)
                    acid_temp = copy.deepcopy(acid)
                    acid_temp.check_interactions(protein_temp)
                    score_direction = protein_temp.score

                    if score_direction < best_score:
                        print("inside")
                        best_score = score_direction
                        best_direction = direction
                        #best_fold_direction = copy.deepcopy(protein)
                    print('location_dir', acid.location)
                    #protein.used_coordinates.remove((tuple(acid.location)))

            #protein.create_bond(acid, protein.sequence_list[i - 1], best_direction)
            
            protein.create_bond(acid, protein.sequence_list[i - 1], best_direction)       
            best_directions.append(best_direction)
            acid.check_interactions(protein)
            print('location', acid.location)
            coordinates = tuple(acid.location)
            used_coordinates_greedy.add(coordinates)
            #protein.used_coordinates.add((tuple(acid.location)))
        # protein.score = 0
        # for acid in protein.sequence_list:
        #     acid.check_interactions(protein)

        print('used_cor', protein.used_coordinates)
    return best_directions


