import math

class Protein():
    """
    A class to represent a protein.
    ...
    Attributes:
    -----------
    sequence: string
        sequence of aminoacids
    sequence_list: list
        list of all aminoacids in this protein
    score: int
        score representing the stability of the protein
    hh_bonds: list of tuples
        contains pairs of aminoacids that have HH-bonds
    ch_bonds: list of tuples
        contains pairs of aminoacids that have CH-bonds
    cc_bonds: list of tuples
        contains pairs of aminoacids that have CC-bonds

    Methods:
    -----------
    add_aminoacid():
        add aminoacid to list
    create_bonds():
        create bond between two aminoacids
    create_output():
        create asked output
    """
    def __init__(self, sequence):
        self.sequence = sequence
        self.sequence_list = []

        self.score = 0
        self.used_coordinates = set()
        
        self.hh_bonds = []
        self.ch_bonds = []
        self.cc_bonds = []
        
    def add_aminoacid(self, acid):
        """
        Add an aminoacid to the protein.
        """
        acid = acid.lower()

        # add the appropriate type of aminoacid based on the sequence
        if acid == 'p':
            self.sequence_list.append(Aminoacid('P'))

        elif acid == 'h':
            self.sequence_list.append(Aminoacid('H'))

        elif acid == 'c':
            self.sequence_list.append(Aminoacid('C'))

        else:
            print("SEQUENCE ERROR: Please only insert aminoacids of type 'H', 'P' or 'C'")

        # initialize the previous neighbour of aminoacid from the second aminoacid in the protein
        if len(self.sequence_list) > 1:
            self.sequence_list[-1].neighbour1 = self.sequence_list[-2]

    def create_bond(self, acid, previous_acid, direction):
        """
        For one aminoacid, determine the direction of the bond, and check for possible interactions.
        """
        # change of direction on x and y axis and step of previous acid
        direction_x = direction[0]
        direction_y = direction[1]
        previous_acid.step = direction[2]

        # create bond between previous and new acid based on step of previous acid
        acid.location_x = previous_acid.location_x + direction_x
        acid.location_y = previous_acid.location_y + direction_y

        coordinates = (acid.location_x, acid.location_y)

        # only validate location if there is no other acid on that location
        if coordinates not in self.used_coordinates:
            acid.location_valid = True
            
        self.used_coordinates.add(coordinates)

    def create_output(self):
        """
        Create output in the asked format.
        """
        print('amino,fold')

        for acid in self.sequence_list:
            print(f'{acid.type},{acid.step}')

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
    type: string
        polar, hydrophobic or cysteine aminoacid
    color: string
        polar = 'royalgreen'
        hydrophobic = 'red'
        cysteine = 'limegreen'

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
        self.location = None
        self.neighbour1 = None

        self.step = 0
        self.type = type

        self.location_valid = False

        # initialize color based on type
        if self.type == 'P':
            self.color = 'b'
        
        elif self.type == 'H':
            self.color = 'r'
        
        elif self.type == 'C':
            self.color = 'g'
    
    def distance(self, other, 3D = False):
        """
        Calculate euclidian distance between this aminoacid and another aminoacid.
        """
        if 3D = True:
            return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        else:
            return math.sqrt((self.location_x - other.location_x) ** 2 + (self.location_y - other.location_y) ** 2)

    def check_interactions(self, protein):
        """
        Check surrounding of aminoacid for other aminoacids, if they are present, check type
        and change score according to the interaction type.
        """
        # loop through all potential interactors for this aminoacid
        for potential_interactor in protein.sequence_list:

            # only check for interactions with other acid than bound amino acids
            if potential_interactor != self.neighbour1:
                
                # only pick interactors in close proximity of acid
                if self.distance(potential_interactor) == 1:
                    location_acid = (self.location_x, self.location_y)
                    location_interactor = (potential_interactor.location_x, potential_interactor.location_y)

                    # change score according to which interaction is happening
                    if self.type == 'H' and potential_interactor.type == 'H':
                        protein.score -= 1
                        protein.hh_bonds.append((location_acid, location_interactor))
                
                    elif self.type == 'H' and potential_interactor.type == 'C':
                        protein.score -= 1
                        protein.ch_bonds.append((location_acid, location_interactor))
                    
                    elif self.type == 'C' and potential_interactor.type == 'H':
                        protein.score -= 1
                        protein.ch_bonds.append((location_acid, location_interactor))
             
                    elif self.type == 'C' and potential_interactor.type == 'C':
                        protein.score -= 5
                        protein.cc_bonds.append((location_acid, location_interactor))