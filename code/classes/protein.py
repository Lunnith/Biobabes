import matplotlib.pyplot as plt
import random
import math

class Protein():
    """
    A class to represent a protein.
    ...
    Attributes:
    -----------
    sequence_list: list
        list of all aminoacids in this protein
    score: int
        score representing the stability
    hh_bonds: list of tuples
        contains pairs of aminoacids that have HH-bonds
    ch_bonds: list of tuples
        contains pairs of aminoacids that have CH-bonds
    cc_bonds: list of tuples
        contains pairs of aminoacids that have CC-bonds

    Methods:
    -----------
    add_aminoacid():
        add aminoacids to list
    initialize_neighbours():
        assign neighbours to aminoacida
    create_bonds():
        create bonds between aminoacids and check possible interactions
    create_output():
        create asked output
    """
    def __init__(self, sequence):
        self.sequence_list = []
        self.score = 0
        self.hh_bonds = []
        self.ch_bonds = []
        self.cc_bonds = []

    def add_aminoacid(self, acid):
        """
        Filling the list of all aminoacids, using the string given on initiation.
        """
        acid = acid.lower()

        if acid == 'p':
            self.sequence_list.append(Aminoacid('P'))

        elif acid == 'h':
            self.sequence_list.append(Aminoacid('H'))

        elif acid == 'c':
            self.sequence_list.append(Aminoacid('C'))
    
    def initialize_neighbours(self):
        """
        Assign neighbours to all aminoacids based on sequence.
        """
        for i in range(len(self.sequence_list)):
            self.sequence_list[i].neighbour1 = self.sequence_list[i - 1]

            if i < len(self.sequence_list) - 1:
                self.sequence_list[i].neighbour2 = self.sequence_list[i + 1]

        self.sequence_list[0].neighbour1 = None

    def create_bonds(self):
        """
        For each aminoacid, determine the direction of the bond, 
        and check for possible interactions.
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
                    else: 
                        amount_of_tries += 1
                    
                used_coordinates.add(coordinates)

                acid.check_interactions(self)
                self.temporary_acids.append(acid)

                previous_acid = acid

    def create_output(self):
        """
        Create output in the asked format.
        """
        print('amino,fold')

        for aminoacid in self.sequence_list:
            print(f'{aminoacid.type},{aminoacid.step}')

        print(f'score,{self.score}')


class Aminoacid():
    """
    A class to represent an aminoacid.
    ...
    Attributes:
    -----------
    location_x: int
        the x coordinate of the position of the aminoacid after folding
    location_y: int
        the y coordinate of the position of the aminoacid after folding
    step: int
        the direction of the folding for this acid

    Methods:
    -----------
    distance():
        calculate distance between one aminoacid and another
    check_interactions():
        check for possible other interacting aminoacids
    """
    def __init__(self, type):
        self.location_x = None
        self.location_y = None
        self.step = 0
        self.type = type

        self.neighbour1 = None
        self.neighbour2 = None

        if self.type == 'P':
            self.color = 'royalgreen'
        
        elif self.type == 'H':
            self.color = 'red'
        
        elif self.type == 'C':
            self.color = 'limegreen'
    

    def distance(self, other):
        """
        Calculate distance between this aminoacid and another aminoacid.
        """
        return math.sqrt((self.location_x - other.location_x) ** 2 + (self.location_y - other.location_y) ** 2)

    def check_interactions(self, protein):
        """
        Check surrounding of aminoacid for other aminoacids, if they are present, check type
        and change score according to the interaction type.
        """
        for acid in protein.temporary_acids:
            if acid != self.neighbour1 and acid != self.neighbour2:
 
                if self.distance(acid) == 1:
                    potential_interactor = acid
                    location_acid = (self.location_x, self.location_y)
                    location_interactor = (acid.location_x, acid.location_y)

                    if self.type == 'H' and potential_interactor.type == 'H':
                        protein.score -= 1
                        protein.hh_bonds.append((location_acid, location_interactor))
                
            
                    if self.type == 'H' and potential_interactor.type == 'C':
                        protein.score -= 1
                        protein.ch_bonds.append((location_acid, location_interactor))
                 
                    
                    if self.type == 'C' and potential_interactor.type == 'H':
                        protein.score -= 1
                        protein.ch_bonds.append((location_acid, location_interactor))
             
            
                    if self.type == 'C' and potential_interactor.type == 'C':
                        protein.score -= 5
                        protein.cc_bonds.append((location_acid, location_interactor))
              
            potential_interactor = None
    