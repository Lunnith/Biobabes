import copy
from ..classes.protein import Protein

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

        self.number_of_states_checked = []
        self.number_of_states = 0
        self.best_scores_list = []

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
                self.number_of_states += 1

    def check_folding(self, new_fold):
        """
        Checks and accepts foldings with higher stability and a lower score than the current best folding.
        """
        new_score = new_fold.score
        old_score = self.best_score

        if new_score <= old_score:
            self.best_score = new_score
            self.best_state = new_fold

            self.best_scores_list.append(self.best_score)
            self.number_of_states_checked.append(self.number_of_states)
    
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
                if directions_pruning and len(new_fold.sequence) >= 5:
                    if len(new_fold.sequence_list) >= (len(new_fold.sequence) * (2/3)) and self.number_of_used_directions(new_fold) < 3:
                        continue

                self.create_child_states(new_fold, new_fold.sequence_list[-1], P_pruning)
            
            else:
                self.check_folding(new_fold)
            
            self.protein = self.best_state   
    
class RecognizePatterns(DepthFirst):
    """
    """
    def recognize_patterns(self):
        """
        """
        P_count = 0
        sequence_parts = []
        sequence_part = []

        for i in range(len(self.protein.sequence)):
            sequence_part.append(self.protein.sequence[i])

            if self.protein.sequence[i] == 'P':
                P_count += 1

            if P_count == 2:
                sequence_parts.append(sequence_part)

                if i != len(self.protein.sequence) - 1:
                    sequence_part = []
                    P_count = 0

        sequence_parts.append(sequence_part)
        return sequence_parts
    
    def depth_fold_in_parts(self, sequence_parts):
        """
        """
        directions = []

        for part in sequence_parts:
            protein_part = Protein(part)
            depth_first = DepthFirst(protein_part, 3)
            print(part)
            depth_first.run(directions_pruning=True)
            print(depth_first.protein)

            steps = []

            for acid in depth_first.protein.sequence_list:
                for direction in self.directions:
                    if direction[3] == acid.step:
                        steps.append(direction)
            
            directions.append(steps)

        return directions

    def connect_folded_parts(self, directions, connection_direction):
        """
        """
        new_protein = Protein(self.protein.sequence)
        all_directions = []

        for directions_part in directions:
            all_directions.extend(directions_part)
            all_directions.append(connection_direction)

        for i in range(len(new_protein.sequence)):
            new_protein.add_aminoacid(new_protein.sequence[i])
            acid = new_protein.sequence_list[i]
        
            # location of first aminoacid is (0,0,0)
            if i == 0:
                acid.location = [0,0,0]
                new_protein.used_coordinates.add((tuple(acid.location)))

            # for other aminoacids then the first create bond with previous acid
            elif i != 0 or i != len(new_protein.sequence) - 1:
                new_protein.create_bond(acid, new_protein.sequence_list[i - 1], all_directions[i - 1])

            # only start checking interactions after the third aminoacid has been added
            if i > 2:
                acid.check_interactions(new_protein)
            
            self.protein = new_protein

    def try_connection_directions(self, directions):
        best_score = 1
        best_fold = None

        for connection_direction in self.directions:
            self.connect_folded_parts(directions, connection_direction)
            score = self.protein.score

            if score < best_score:
                best_fold = self.protein

        return best_fold

    def run_in_parts(self):
        """
        """
        sequence_parts = self.recognize_patterns()
        directions = self.depth_fold_in_parts(sequence_parts)
        best_fold = self.try_connection_directions(directions)

        self.protein = best_fold


# 14.55 start
# 36 aminozuren

# to do: location_valid (niet overzichzelf heen en vastlopen) implementeren
# to do: voor elke connection los elke richting kunnen proberen
        









