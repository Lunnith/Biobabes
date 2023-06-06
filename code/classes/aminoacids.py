import random

class Protein():
    """
    
    """
    def __init__(self, sequence):
        """
        Initiate the dictionary of all aminoacids.
        For structure of this dictionary, see add_aminoacid.

        Initiate the possible directions.
        """
        self.sequence_dict = {}
        self.add_aminoacid(sequence)

        self.directions = [-1, 1, 2, -2]

    def add_aminoacid(self, sequence):
        """
        Filling the dictionary of all aminoacids, using the string given on initiation.

        The keys of this dictionary are the created aminoacid objects, 
        The values are set on 'None', for now, but will later on be set to the direction of the bond.
        """
        for acid in sequence.lower():
            if acid == 'p':
                self.sequence_dict[Polar()] = None
            elif acid == 'h':
                self.sequence_dict[Hydrophobic()] = None
            elif acid == 'c':
                self.sequence_dict[Cysteine()] = None
    
    def create_bonds(self):
        """
        For each aminoacid, determine the direction of the bond.
        """
        for acid in self.sequence_dict:
            self.sequence_dict[acid] = random.choice(self.directions)

    
    def check_bonds(self):
        """
        """
        pass

    def create_output(self):
        """
        """
        pass




class Aminoacid():
    """
    
    """
    def __init__(self):
        """
        """
        pass

    def interact(self, other):
        """
        """
        pass
        
    
class Polar(Aminoacid):
    """
    Create a polar aminoacid.
    """
    def __init__(self):
        """
        """
        super().__init__()

        self.color = 'royalblue'

class Hydrophobic(Aminoacid):
    """
    Create a Hydrophobic aminoacid.
    """
    def __init__(self):
        super().__init__()

        self.color = 'red'

class Cysteine(Aminoacid):
    """
    Create a Cysteine aminoacid.
    """
    def __init__(self):
        super().__init__()

        self.color = 'lime'
    


protein = Protein("HPCHP")
protein.create_bonds()
print(protein.sequence_dict)