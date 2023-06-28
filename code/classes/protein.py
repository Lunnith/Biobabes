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
    add_aminoacid(acid):
        add aminoacid to list
    create_bond(acid, previous_acid, direction):
        create bond between two aminoacids
    get_temp_sequence(temp_protein):
        returns temporary sequence of temporary protein
    is_valid():
        checks if all locations of acids are valid
    create_output():
        create asked output
    """
    def __init__(self, sequence: str) -> None:
        allowed_types = set(('H', 'P', 'C'))
        sequence_set = set(sequence)

        for type in sequence_set:
            if type not in allowed_types:
                raise ValueError("SEQUENCE ERROR: Please only insert aminoacids of type 'H', 'P' or 'C'")
        
        self.sequence = sequence
        self.sequence_list = []

        self.score = 0
        self.used_coordinates = set()
        
        self.hh_ch_bonds = []
        self.cc_bonds = []

        self.valid = True
        
    def add_aminoacid(self, acid: str) -> None:
        """
        Add an aminoacid to the protein.
        """
        self.sequence_list.append(Aminoacid(acid))

        # initialize the previous neighbour of aminoacid from the second aminoacid in the protein
        if len(self.sequence_list) > 1:
            self.sequence_list[-1].neighbour1 = self.sequence_list[-2]

    def create_bond(self, acid: Aminoacid, previous_acid: Aminoacid, direction: tuple) -> None:
        """
        For one aminoacid, determine the direction of the bond..
        """
        previous_acid.step = direction[3]

        # create bond between previous and new acid based on step of previous acid
        acid.location = list(map(add, previous_acid.location, direction[0:3]))

        coordinates = tuple(acid.location)

        # only validate location if there is no other acid on that location
        if coordinates not in self.used_coordinates:
            acid.location_valid = True
            
            self.used_coordinates.add(coordinates)

    def get_temp_sequence(self) -> str:
        """
        Method to get the temporary sequence of a protein that is being build.
        """
        temp_sequence = str()

        for acid in self.sequence_list:
            temp_sequence += acid.type
        
        return temp_sequence
    
    def is_valid(self) -> bool:
        """
        Method to check if protein didn't fold over itself.
        """
        for acid in self.sequence_list:

            if acid.location_valid == False:
                self.valid = False
        
        return self.valid
    
    def create_output(self) -> None:
        """
        Method to create output in the asked format.
        """
        print('amino,fold')

        for acid in self.sequence_list:
            print(f'{acid.type},{acid.step}')

        print(f'score,{self.score}')