import copy

class DepthFirst():
    """
    A Depth First algorithm that builds a stack of the same protein with a unique folding for each instance.
    """
    def __init__(self, protein, dimensions):
        """
        Method to initialize attributes and methods of class DepthFirst.
        """
        # initialize protein and add first aminoacid with location and step
        self.protein = protein
        self.protein.add_aminoacid(self.protein.sequence[0])
        self.protein.sequence_list[0].location = [0,0,0]
        self.protein.used_coordinates.add((tuple([0,0,0])))

        # initialize directions based on dimensions
        if dimensions == 2:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
        if dimensions == 3:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))

        self.states = [copy.deepcopy(self.protein)]
        self.best_state = None
        self.best_score = 0

    def get_next_state(self):
        """
        Method that gets the next state from the list of states.
        """
        return self.states.pop()
    
    def create_child_states(self, temp_protein, last_added_aminoacid, P_pruning=False):
        """
        Creates all possible child-states with different folding directions and adds them to the list of states. 
        Ignores illegal states. Optional non-optimal pruning that ensures that several P's in a row don't fold in the same direction.
        """
        for direction in self.directions:
            new_fold = copy.deepcopy(temp_protein)
            to_add_aminoacid_type = new_fold.sequence[len(new_fold.get_temp_sequence())]

            # expand new fold with aminoacid and create a bond in the given direction
            new_fold.add_aminoacid(to_add_aminoacid_type)
            new_fold.create_bond(new_fold.sequence_list[-1], new_fold.sequence_list[-2], direction)

            # non-optimal pruning which makes sure that several P's in a row don't fold in the same direction
            if P_pruning and len(new_fold.sequence_list) > 2:
                if new_fold.sequence_list[-1].type == 'P' and new_fold.sequence_list[-2].type == 'P' and new_fold.sequence_list[-2].step == new_fold.sequence_list[-3].step:
                    continue

            # only add new folding to the list if the folding is valid
            if new_fold.sequence_list[-1].location_valid == True:
                new_fold.sequence_list[-1].check_interactions(new_fold)
                self.states.append(new_fold)

    def check_folding(self, new_fold):
        """
        Checks and accepts foldings with higher stability and a lower score than the current best folding.
        """
        new_score = new_fold.score
        old_score = self.best_score

        if new_score <= old_score:
            self.best_score = new_score
            self.best_state = new_fold
    
    def number_of_used_directions(self, temp_protein):
        """
        Method that counts the number of used directions during folding of the protein.
        """
        used_directions = set()

        for acid in temp_protein.sequence_list:
            used_directions.add(acid.step)

        return len(used_directions)

    def run(self, P_pruning=False, directions_pruning=False):
        """
        Runs the Depth First Algorithm until it has seen all possible states.
        """
        depth = len(self.protein.sequence)

        # loop through stack of foldings
        while self.states:
            new_fold = self.get_next_state()

            # if protein is not complete yet, expand the folding, else check the score of the new folding
            if len(new_fold.sequence_list) < depth:
                
                # non-optimal pruning that skips states that have used less than 3 directions at 2/3 of the folding
                if directions_pruning:
                    if len(new_fold.sequence_list) >= (len(new_fold.sequence) * (2/3)) and self.number_of_used_directions(new_fold) < 3:
                        continue

                self.create_child_states(new_fold, new_fold.sequence_list[-1], P_pruning)
            
            else:
                self.check_folding(new_fold)
            
            self.protein = self.best_state   
    
    



