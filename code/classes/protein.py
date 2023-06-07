import matplotlib.pyplot as plt
import random
import math

class Protein():
    """
    
    """
    def __init__(self, sequence):
        """
        Initiate the dictionary of all aminoacids.
        For structure of this dictionary, see add_aminoacid.

        Initiate the possible directions.
        """
        self.sequence_list = []

        self.add_aminoacid(sequence)
        self.initialize_neighbours()
        self.score = 0


    def add_aminoacid(self, sequence):
        """
        Filling the dictionary of all aminoacids, using the string given on initiation.

        The keys of this dictionary are the created aminoacid objects, 
        The values are set on 'None', for now, but will later on be set to the direction of the bond.
        """
        for acid in sequence.lower():
            if acid == 'p':
                self.sequence_list.append(Polar())
            elif acid == 'h':
                self.sequence_list.append(Hydrophobic())
            elif acid == 'c':
                self.sequence_list.append(Cysteine())
    
    def initialize_neighbours(self):
        """
        
        """
        for i in range(len(self.sequence_list)):
            self.sequence_list[i].neighbour1 = self.sequence_list[i - 1]
            if i < len(self.sequence_list) - 1:
                self.sequence_list[i].neighbour2 = self.sequence_list[i + 1]

        self.sequence_list[0].neighbour1 = None

    def create_bonds(self):
        """
        For each aminoacid, determine the direction of the bond.
        """
        self.sequence_list[0].location_x = 0
        self.sequence_list[0].location_y = 0
        previous_acid = self.sequence_list[0]
        directions = [[1, 0, 1], [-1, 0, -1], [0, 1, 2], [0, -1, -2]]

        self.temporary_acids = []
        used_coordinates = set()

        for acid in self.sequence_list:
            
            if acid == self.sequence_list[0]:
                self.temporary_acids.append(acid)
                used_coordinates.add((acid.location_x, acid.location_y))
                

            else:
                location_valid = False
                amount_of_tries = 0

                while location_valid == False:

                    # Select a direction with change in x and y cordinates
                    direction = random.choice(directions)
                    
                    # define direction that x and y will go in
                    direction_x = direction[0]
                    direction_y = direction[1]

                    # Define step for the acid
                    previous_acid.step = direction[2]

                    # Set location of the acid
                    acid.location_x = previous_acid.location_x + direction_x
                    acid.location_y = previous_acid.location_y + direction_y

                    coordinates = (acid.location_x, acid.location_y)

                    if coordinates not in used_coordinates and amount_of_tries < 20:
                        location_valid = True
                    else: amount_of_tries += 1
                    
                used_coordinates.add(coordinates)

                acid.check_interactions(self)
                self.temporary_acids.append(acid)

                previous_acid = acid

    def create_output(self):
        """
        """
        print('amino,fold')
        for aminoacid in self.sequence_list:
            print(f'{aminoacid.type},{aminoacid.step}')
        print(f'score,{self.score}')


class Aminoacid():
    """
    
    """
    def __init__(self):
        """
        """
        self.location_x = None
        self.location_y = None
        self.step = 0

        
        self.neighbour1 = None
        self.neighbour2 = None

    def distance(self, other):
        return math.sqrt((self.location_x - other.location_x) ** 2 + (self.location_y - other.location_y) ** 2)

    def check_interactions(self, protein):
        """
        Check surrounding of aminoacid for other aminoacids, if they are present, check type
        and change score according to the interaction type.
        """

        for acid in protein.temporary_acids:
            if acid != self.neighbour1 or acid != self.neighbour2:
                print(self.distance(acid))
                if math.isclose(self.distance(acid), 1):
                    print('if', self.distance(acid))
                    potential_interactor = acid
                    print('pot', potential_interactor)
                    print('t self', type(self))
                    print('t pot', type(potential_interactor))
                    if type(self) == Hydrophobic and type(potential_interactor) == Hydrophobic:
                        protein.score -= 1
                        print("erin")
            
                    if type(self) == Hydrophobic and type(potential_interactor) == Cysteine:
                        protein.score -= 1
                        print("erin")
                    
                    if type(self) == Cysteine and type(potential_interactor) == Hydrophobic:
                        protein.score -= 1
                        print("erin")
            
                    if type(self) == Cysteine and type(potential_interactor) == Cysteine:
                        protein.score -= 5
                        print("erin")
            
            potential_interactor = None
        
    
class Polar(Aminoacid):
    """
    Create a polar aminoacid.
    """
    def __init__(self):
        """
        """
        super().__init__()

        self.color = 'royalblue'
        self.type = 'P'
    

class Hydrophobic(Aminoacid):
    """
    Create a Hydrophobic aminoacid.
    """
    def __init__(self):
        super().__init__()

        self.color = 'red'
        self.type = 'H'
     

class Cysteine(Aminoacid):
    """
    Create a Cysteine aminoacid.
    """
    def __init__(self):
        super().__init__()

        self.color = 'limegreen'
        self.type = 'C'
     
    