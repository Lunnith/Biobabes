import random
import copy

def random_assignment(protein, directions):
    """
    Randomly assign direction to each bond within the protein.
    """
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
