from .aminoacid import Aminoacid
from operator import add

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
    hh_ch_bonds: list of tuples
        contains pairs of aminoacids that have HH-bonds or CH-bonds
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
        allowed_types = set(('H', 'P', 'C'))
        sequence_set = set(sequence)

        if allowed_types != sequence_set:
            raise ValueError("SEQUENCE ERROR: Please only insert aminoacids of type 'H', 'P' or 'C'")
        
        self.sequence = sequence
        self.sequence_list = []

        self.score = 0
        self.used_coordinates = set()
        
        self.hh_ch_bonds = []
        self.cc_bonds = []
        
    def add_aminoacid(self, acid):
        """
        Add an aminoacid to the protein.
        """
        self.sequence_list.append(Aminoacid(acid))

        # initialize the previous neighbour of aminoacid from the second aminoacid in the protein
        if len(self.sequence_list) > 1:
            self.sequence_list[-1].neighbour1 = self.sequence_list[-2]

    def create_bond(self, acid, previous_acid, direction):
        """
        For one aminoacid, determine the direction of the bond..
        """
        previous_acid.step = direction[3]

        # create bond between previous and new acid based on step of previous acid
        acid.location = list(map(add, previous_acid.location, direction[0:3]))

        #coordinates = (acid.location_x, acid.location_y, acid.location_z)
        coordinates = tuple(acid.location)

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