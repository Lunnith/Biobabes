from .depth_first import DepthFirst
from ..classes.protein import Protein

import random
import copy

class ImportantParts(DepthFirst):
    """
    A class to separate a big protein into smaller parts, which applies the depth first algorithm to all separate parts.
    ...
    Attributes:
    -----------
    Inherits all attributes from DepthFirst.

    Methods:
    -----------
    Inherits all methods from DepthFirst. Extra methods:

    recognize_important_parts(split_on_P, split_on_size, size):
        splits the protein into parts, based on P count or size
    depth_first_in_parts(sequence_parts):
        performs depth first for alle separate parts
    connect_folded_parts(best_directions):
        connects all folded parts
    run(iterations, split_on_P, split_on_size, size):
        runs important part algorithm a number of times and saves the best folding

    """
    def recognize_important_parts(self, split_on_P: bool, split_on_size: bool, size: int) -> list:
        """
        Method to recognize important parts in the sequence and separate the protein into smaller parts, based on the
        count of P's or the length.
        """
        if split_on_P and split_on_size:
            raise Exception("You can't split on size and on P_count")
        
        P_count = 0
        sequence_part = []
        sequence_parts_list = []
        
        # loop through all aminoacids in protein sequence and separate sequence into smaller parts based on P count or size
        for i in range(len(self.protein.sequence)):
            sequence_part.append(self.protein.sequence[i])

            if split_on_P:
                if self.protein.sequence[i] == 'P':
                    P_count += 1

                # cut protein if P count is 2 or if part is bigger than 10 aminoacids
                if P_count == 2 or len(sequence_part) > 10:
                    sequence_parts_list.append(sequence_part)

                    if i != len(self.protein.sequence) - 1:
                        sequence_part = []
                        P_count = 0
            
            if split_on_size:
                if len(sequence_part) == size:
                    sequence_parts_list.append(sequence_part)

                    if i != len(self.protein.sequence) - 1:
                        sequence_part = []

        sequence_parts_list.append(sequence_part)
        return sequence_parts_list
    
    def depth_first_in_parts(self, sequence_parts: list) -> list:
        """
        Method to apply depth first algorithm on all separate parts.
        """
        best_directions = []

        # loop through sequence parts and apply depth first to every part
        for part in sequence_parts:
            protein_part = Protein(part)
            depth_first = DepthFirst(protein_part, self.dimensions)
            depth_first.run(directions_pruning=True, P_pruning=True)

            part_directions = []

            # save the directions for the optimal folding of every part
            for acid in depth_first.protein.sequence_list:
                for direction in self.directions:
                    if direction[3] == acid.step:
                        part_directions.append(direction)

            best_directions.append(part_directions)

        return best_directions

    def connect_folded_parts(self, best_directions: list) -> Protein:
        """
        Method to connect the optimal folding of all parts of the protein.
        """
        protein = Protein(self.protein.sequence)
        all_directions = []

        # create list of all best directions and connection directions between the different parts for the final protein
        for directions_of_part in best_directions:
            all_directions.extend(directions_of_part)
            all_directions.append(random.choice(list(self.directions)))

        new_protein = copy.deepcopy(protein)

        # create protein based on best directions
        for i in range(len(new_protein.sequence)):
            new_protein.add_aminoacid(new_protein.sequence[i])
            acid = new_protein.sequence_list[i]

            # first aminoacid has a defined location
            if i == 0:
                acid.location = [0,0,0]
                new_protein.used_coordinates.add((tuple(acid.location)))
                acid.location_valid = True

            # for other aminoacids than the first and the last, create bond with previous acid and check interactions
            elif i != 0:
                new_protein.create_bond(acid, new_protein.sequence_list[i - 1], all_directions[i - 1])
                acid.check_interactions(new_protein)
            
        return new_protein

    def run(self, iterations: int = 1000, split_on_P: bool = True, split_on_size: bool = False, size: int = 12) -> None:
        """
        Method to run the depth first algorithm for separate parts of the protein and connect these parts
        """
        # split the protein and make a list of the directions of all optimal foldings
        sequence_parts = self.recognize_important_parts(split_on_P, split_on_size, size)
        best_directions = self.depth_first_in_parts(sequence_parts)

        best_fold = None
        best_score = 1

        # find the protein with the highest score based on the different random connection directions that is valid
        for i in range(iterations):
            new_protein = self.connect_folded_parts(best_directions)
            score = new_protein.score

            if score < best_score and new_protein.is_valid():
                best_fold = new_protein

        if best_fold == None:
            self.protein = new_protein
            self.protein.score = 0
        else:
            self.protein = best_fold