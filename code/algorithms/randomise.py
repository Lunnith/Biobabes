import random
import copy
from operator import add

def random_assignment(protein, dimensions):
    """
    Randomly assign direction to each bond within the protein.
    """
    if dimensions == 2:
        directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))

    if dimensions == 3:
        directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))

    # loop through protein sequence and add aminoacids accordingly
    for i in range(len(protein.sequence)):
        protein.add_aminoacid(protein.sequence[i])
        acid = protein.sequence_list[i]
        
        # location of first aminoacid is (0,0,0)
        if i == 0:
            acid.location = [0,0,0]
        
            protein.used_coordinates.add((tuple(acid.location)))

        # for other aminoacids then the first create bond with previous acid
        elif i != 0:
            tried_directions = set()

            while acid.location_valid == False:
                direction = random.choice(tuple(directions))
                protein.create_bond(acid, protein.sequence_list[i - 1], direction)
                tried_directions.add(direction)

                # if protein can't fold anymore, return shorter folded protein
                if tried_directions == directions:
                    return protein

        # only start checking interactions after the third aminoacid has been added
        if i > 2:
            acid.check_interactions(protein)

    return protein

def random_reassignment(protein, dimensions, k=20):
    """
    Algorithm that randomly tries several configurations of the same protein and saves the best
    configuration with the highest stability.
    """
    best_score =  1
    scores = []

    # randomly fold the protein k times and save the ones that are better than the previously saved
    for i in range(k):
        protein.__init__(protein.sequence)
        folded_protein = random_assignment(protein, dimensions)

        # only check the score of proteins that were folded completely
        if len(folded_protein.sequence_list) == len(protein.sequence):        
            current_score = folded_protein.score
            scores.append(current_score)

            if current_score < best_score:
                best_score = current_score
                best_fold = copy.deepcopy(folded_protein)
    
    return best_fold, scores
