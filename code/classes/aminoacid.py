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

        self.color = 'lime'
        self.type = 'C'
    