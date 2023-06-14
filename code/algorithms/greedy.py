# from code.classes.protein import Protein
# from code.visualization.visualize import *
import random
import copy
import matplotlib.pyplot as plt
from operator import add

# class Greedy():

#     def __init__(self, splits):

#         self.splits = splits

def greedy(protein):
    
    best_score = 0
    best_directions = [] 
   # if dimensions == 2:
    directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))

    # if dimensions == 3:
    #     directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))

    for i in range(len(protein.sequence)):
        protein.add_aminoacid(protein.sequence[i])
        acid = protein.sequence_list[i]
        
        # location of first aminoacid is (0,0,0)
        if i == 0:
            acid.location = [0, 0, 0]
        
            protein.used_coordinates.add((tuple(acid.location)))

        # for other aminoacids then the first create bond with previous acid
        elif i != 0:
          #  tried_directions = set()
            best_direction = random.choice(tuple(directions))
            for direction in directions:
                protein.score = 0
                protein.create_bond(acid, protein.sequence_list[i - 1], direction)
                if acid.location_valid == False:
                    pass
                else:
                    acid.check_interactions(protein)
                    score_direction = protein.score

                    if score_direction < best_score:
                        print("inside")
                        best_score = score_direction
                        best_direction = direction
                        best_fold_direction = copy.deepcopy(protein)
            
            tried_directions = set()
            while acid.location_valid == False:

                best_direction = random.choice(tuple(directions))
                protein.create_bond(acid, protein.sequence_list[i - 1], direction)
                tried_directions.add(best_direction)

                # if protein can't fold anymore, return shorter folded protein
                if tried_directions == directions:
                    return protein
                
            best_directions.append(best_direction)
            protein.create_bond(acid, protein.sequence_list[i - 1], best_direction)
            acid.check_intercation(protein)
    return best_directions


# proteinD = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
# best_directions = greedy(proteinD, 2)
# print(best_directions)
# proteinD.create_output()
# visualize_protein(proteinD, 2)