from __future__ import annotations
import math

class Aminoacid():
    """
    A class to represent an aminoacid.
    ...
    Attributes:
    -----------
    location: list
        the coordinates of the position of the aminoacid after folding
    step: int
        the direction of the folding for this acid
    type: string
        polar, hydrophobic or cysteine aminoacid
    color: string
        polar = 'blue'
        hydrophobic = 'red'
        cysteine = 'green'

    Methods:
    -----------
    distance(other):
        calculate distance between one aminoacid and another
    check_interactions(protein):
        check for possible other interacting aminoacids
    """
    def __init__(self, type: str) -> None:
     
        self.location = [None, None, None]
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
    
    def distance(self, other: Aminoacid) -> float:
        """
        Calculate euclidian distance between this aminoacid and another aminoacid.
        """
        return math.sqrt((self.location[0] - other.location[0])**2 + (self.location[1] - other.location[1])**2 + (self.location[2] - other.location[2])**2)
   
    def check_interactions(self, protein: Protein, index: int = -1) -> None:
        """
        Check surrounding of aminoacid for other aminoacids, if they are present, check type
        and change score according to the interaction type.
        """
        # loop through all potential interactors for this aminoacid
        for potential_interactor in protein.sequence_list[:index]:

            # only check for interactions with other acid than bound amino acids and interactors in close proximity
            if potential_interactor != self.neighbour1 and self.distance(potential_interactor) == 1:
                
                # change score according to which interaction is happening
                types = set((('H', 'H'), ('C', 'H'), ('H', 'C'), ('C', 'C')))
                interaction_combination = tuple((self.type, potential_interactor.type))

                if interaction_combination in types:
                    
                    if interaction_combination == ('C', 'C'):
                        protein.score -= 5
                        protein.cc_bonds.append((self.location, potential_interactor.location))
                    else:
                        protein.score -= 1
                        protein.hh_ch_bonds.append((self.location, potential_interactor.location))