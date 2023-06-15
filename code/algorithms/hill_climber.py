from .randomise import random_assignment
import random

class Hill_climber():
    def __init__(self, protein, dimensions=3):
        self.protein = protein

        #Initiate first random folding
        self.protein = random_assignment(protein, dimensions)
    
    def change_one_bond(self, protein, dimensions=3):
        #Set options for dimensions
        if dimensions == 2:
            directions = {1: [1, 0, 0], -1: [-1, 0, 0], 2: [0, 1, 0], -2: [0, -1, 0]}
        if dimensions == 3:
            directions = {1: [1, 0, 0], -1: [-1, 0, 0], 2: [0, 1, 0], -2: [0, -1, 0], 3: [0, 0, 1], -3: [0, 0, -1]}
            
        index_changing_bond = random.randint(1, len(protein.sequence_list)-1)
        location_previous_aminoacid = protein.sequence_list[index_changing_bond - 1].location
        old_location = protein.sequence_list[index_changing_bond].location
        old_step = protein.sequence_list[index_changing_bond].step
        acid = protein.sequence_list[index_changing_bond]

        #Get the old direction
        old_direction = directions[old_step]
        print("Old direction is", old_direction)
        print("Old step is", old_step)

        #Determine new step
        tried_directions = set()
        # tried_directions.add(old_direction)
        acid.location_valid = False
        
        while acid.location_valid == False:
            new_step = random.sample(directions.keys(), 1)[0]
            print("New step is", new_step)
            #Overwrite step of previous acid
            protein.sequence_list[index_changing_bond - 1].step = new_step
            #fix format of directions
            directions[new_step].append(new_step)
            print("New directions of this aminoacid", directions[new_step])
            direction = tuple(directions[new_step])
            print("New direction as tuple:", direction)

            #Create new bond based on new direction
            protein.create_bond(acid, protein.sequence_list[index_changing_bond - 1], direction)
            new_location = protein.sequence_list[index_changing_bond].location
            # tried_directions.add(direction)

            # if protein can't fold anymore, return shorter folded protein
            #This doesnt work for dictionaries like this
            if tried_directions == directions:
                print('ended')
                return protein
                

            print(location_previous_aminoacid, old_location, old_step, new_step, new_location)
    
    # def refold(protein)

