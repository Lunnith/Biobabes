import random
import copy

def random_assignment(protein, directions):
    """
    Randomly assign direction to each bond within the protein.
    """
    folded_protein = copy.deepcopy(protein)

    for i in range(len(folded_protein.sequence)):
        folded_protein.add_aminoacid(folded_protein.sequence[i])

        acid = folded_protein.sequence_list[i]
        
        if i == 0:
            acid.location_x = 0
            acid.location_y = 0
            folded_protein.used_coordinates.add((acid.location_x, acid.location_y))

        if i != 0:
            tried_directions = set()

            while acid.location_valid == False and tried_directions.difference(directions) == set():
                direction = random.choice(tuple(directions))
                folded_protein.create_bond(acid, folded_protein.sequence_list[i - 1], direction)

                tried_directions.add(direction)

    return folded_protein
