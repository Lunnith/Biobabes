from code.classes.protein import Protein
from code.visualization.visualize import *
import random
import matplotlib.pyplot as plt
from operator import add
# proteinA = Protein("HPCHPCHH")
proteinA = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
# print(protein.sequence_list)
# random_assignment(proteinA)

# for acid in proteinA.sequence_list:
#     print(f"Class = {acid}, Step = {acid.step}, Coordinates = ({acid.location_x}, {acid.location_y})")

# proteinA.create_output()
# visualize_protein(proteinA)

listA = [0, 1, 0]
listB = [2, 1, 0]
print(list(map(add, listA, listB)))

def make_and_fold_protein(sequence, k=1):
    protein = Protein(sequence)
    directions = set(((1, 0, 1), (-1, 0, -1), (0, 1, 2), (0, -1, -2)))

    for i in range(len(protein.sequence)):
        protein.add_aminoacid(protein.sequence[i])

        acid = protein.sequence_list[i]
        
        if i == 0:
            acid.location_x = 0
            acid.location_y = 0
            protein.used_coordinates.add((acid.location_x, acid.location_y))

        if i != 0:
            tried_directions = set()

            while acid.location_valid == False and tried_directions.difference(directions) == set():
                direction = random.choice(tuple(directions))
                protein.create_bond(acid, protein.sequence_list[i - 1], direction)

                tried_directions.add(direction)

        if i > 2:
            acid.check_interactions(protein)

    return protein

sequence = 'HPHPCHCPHCH'
protein = make_and_fold_protein(sequence)
protein.create_output()
visualize_protein(protein)


