import random
import copy

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
                    protein.valid = False
                    return protein

        # only start checking interactions after the third aminoacid has been added
        if i > 2:
            acid.check_interactions(protein)

    return protein

def random_reassignment(protein, dimensions, k=20, func = random_assignment):
    """
    Algorithm that randomly tries several configurations of the same protein and saves the best
    configuration with the highest stability.
    """
    best_score =  1
    scores = []

    # randomly fold the protein k times and save the ones that are better than the previously saved
    for i in range(k):
        protein.__init__(protein.sequence)
        folded_protein = func(protein, dimensions)

        # only check the score of proteins that were folded completely
        if folded_protein.valid == True:  
            current_score = folded_protein.score
            scores.append(current_score)

            # check if the score of this protein is better than the best score
            if current_score < best_score:

                # set the current score as the new best score and a copy of the protein as the new best fold
                best_score = current_score
                best_fold = copy.deepcopy(folded_protein)
    
    # return the best fold, all scores in a list and the best score
    return best_fold, scores, best_score
