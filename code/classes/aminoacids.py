import random

class Protein():
    """
    Superclass
    """
    def __init__(self, sequence):
        """
        Initiate the dictionary of all aminoacids, where the keys are 
        """
        self.sequence_dict = {}
        self.add_aminoacid(sequence)

        self.directions = [-1, 1, 2, -2]

    def add_aminoacid(self, sequence):
        """
        """
        for acid in sequence.lower():
            if acid == 'p':
                self.sequence_dict[Polar(acid)] = None
            elif acid == 'h':
                self.sequence_dict[Hydrophobic(acid)] = None
            elif acid == 'c':
                self.sequence_dict[Cysteine(acid)] = None
    
    def create_bonds(self):
        """
        """
        for acid in self.sequence_dict:
            self.sequence_dict[acid] = random.choice(self.directions)

    
    def check_bonds(self):
        """
        """

    def create_output(self):
        """
        """




class Aminoacid():
    """
    
    """
    def __init__(self, type):
        """
        """

    def interact(self, other):
        """
        """
        
    
class Polar(Aminoacid):
    """
    
    """
    def__init__(self):
        """
        """
        super().__init__()

        self.color = 'royalblue'

class Hydrophobic(Aminoacid):
    """
    
    """
    def __init__(self):
        super().__init__()

        self.color = 'red'

class Cysteine(Aminoacid):
    """
    
    """
    def __init__(self, type):
        super().__init__()

        self.color = 'lime'
    


