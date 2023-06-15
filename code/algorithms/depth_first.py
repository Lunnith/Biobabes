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
        self.protein.sequence_list[0].step = 1
        self.protein.used_coordinates.add((tuple([0,0,0])))

        # initialize directions based on dimensions
        if dimensions == 2:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2)))
        else:
            self.directions = set(((1, 0, 0, 1), (-1, 0, 0, -1), (0, 1, 0, 2), (0, -1, 0, -2), (0, 0, 1, 3), (0, 0, -1, -3)))
        
        self.states = [copy.deepcopy(self.protein)]
        self.best_state = None
        self.best_score = 0

    def get_next_state(self):
        """
        Method that gets the next state from the list of states.
        """
        return self.states.pop()
    
    def build_children(self, temp_protein, last_added_aminoacid):
        """
        Creates all possible child-states with different folding directions and adds them to the list of states. 
        Ignores illegal states.
        """
        for direction in self.directions:

            # only make a new state if direction is not the reverse direction of the previous bond
            if direction[3] != -last_added_aminoacid.step:
                new_fold = copy.deepcopy(temp_protein)

                # determine the new to add aminoacid based on the temporary sequence of the old state
                temp_sequence = str()
                for acid in new_fold.sequence_list:
                    temp_sequence += acid.type

                to_add_aminoacid_type = new_fold.sequence[len(temp_sequence)]

                # expand new fold with aminoacid and create a bond in the given direction
                new_fold.add_aminoacid(to_add_aminoacid_type)
                new_fold.create_bond(new_fold.sequence_list[-1], new_fold.sequence_list[-2], direction)

                # only add new folding to the list if the folding is valid
                if new_fold.sequence_list[-1].location_valid == True:
                    new_fold.sequence_list[-1].check_interactions(new_fold)
                    self.states.append(new_fold)

    def check_solution(self, new_fold):
        """
        Checks and accepts better solutions than the current solution.
        """
        new_score = new_fold.score
        old_score = self.best_score

        if new_score < old_score:
            self.best_score = new_score
            self.best_state = new_fold
            print(f'New best score: {self.best_score}')

    def run(self):
        """
        Runs the Depth First Algorithm until it has seen all possible states.
        """
        depth = len(self.protein.sequence)

        # loop through stack of foldings
        while self.states:
            new_fold = self.get_next_state()

            # if protein is not complete yet, expand the folding, else check the score of the new folding
            if len(new_fold.sequence_list) < depth:
                self.build_children(new_fold, new_fold.sequence_list[-1])
            else:
                self.check_solution(new_fold)

            self.protein = self.best_state
