import random
import copy

def random_assignment(protein, directions):
    """
    Randomly assign direction to each bond within the protein.
    """

    # loop through protein sequence and add aminoacids accordingly
    for i in range(len(protein.sequence)):
        protein.add_aminoacid(protein.sequence[i])
        acid = protein.sequence_list[i]
        
        # location of first aminoacid is (0,0)
        if i == 0:
            acid.location_x = 0
            acid.location_y = 0
            protein.used_coordinates.add((acid.location_x, acid.location_y))

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

def random_reassignment(protein, directions, k=20):
    """
    Algorithm that randomly tries several configurations of the same protein and saves the best
    configuration with the highest stability.
    """
    best_score =  1
    scores = []

    # randomly fold the protein k times and save the ones that are better than the previously saved
    for i in range(k):
        protein.__init__(protein.sequence)
        folded_protein = random_assignment(protein, directions)

        # only check the score of proteins that were folded completely
        if len(folded_protein.sequence_list) == len(protein.sequence):        
            current_score = folded_protein.score
            scores.append(current_score)

            if current_score < best_score:
                best_score = current_score
                best_fold = copy.deepcopy(folded_protein)
    
    return best_fold, scores
