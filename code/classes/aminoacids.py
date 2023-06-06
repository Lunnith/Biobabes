import random
import matplotlib.pyplot as plt

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
        self.score = 0
        self.add_aminoacid(sequence)


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
    
    def create_bonds(self):
        """
        For each aminoacid, determine the direction of the bond.
        """
        self.sequence_list[0].location_x = 0
        self.sequence_list[0].location_y = 0
        previous_acid = self.sequence_list[0]
        directions = [[1, 0, 1], [-1, 0, -1], [0, 1, 2], [0, -1, -2]]

        for acid in self.sequence_list:
            
            if acid == self.sequence_list[0]:
                pass

            else:
                # Select a direction with change in x and y cordinates
                direction = random.choice(directions)
                
                # define direction that x and y will go in
                direction_x = direction[0]
                direction_y = direction[1]

                # Define step for the acid
                acid.step = direction[2]

                # Set location of the acid
                acid.location_x = previous_acid.location_x + direction_x
                acid.location_y = previous_acid.location_y + direction_y


            previous_acid = acid


    def check_interactions(self):
        """
        Check surrounding of aminoacid for other aminoacids, if they are present, check type
        and change score according to the interaction type.
        """
        for acid in self.sequence_list:
            if acid.location_x == self.location_x + 1 and acid.location_y == self.location_y:
                other = acid
            elif acid.location_x == self.location_x - 1 and acid.location_y == self.location_y:
                other = acid
            elif acid.location_x == self.location_x and acid.location_y == self.location_y + 1:
                other = acid
            elif acid.location_x == self.location_x and acid.location_y == self.location_y - 1:
                other = acid

        if type(self) == Hydrophobic() and type(other) == Hydrophobic():
            self.score -= 1
        
        if type(self) == Hydrophobic() and type(other) == Cysteine():
            self.score -= 1
        
        if type(self) == Cysteine() and type(other) == Cysteine():
            self.score -= 5

    def create_output(self):
        """
        """
        print('amino,fold')
        for aminoacid in self.sequence_list:
            print(f'{aminoacid.type},{aminoacid.step}')
        print(f'score,{self.score}')
    

    def visualize(self):
        pos_x = []
        pos_y = []
        color = []
        for aminoacid in self.sequence_list:
            pos_x.append(aminoacid.location_x)
            pos_y.append(aminoacid.location_y)
            color.append(aminoacid.color)

        plt.scatter(pos_x, pos_y, c = color)
        plt.plot(pos_x, pos_y, 'b-')
        plt.grid(True)
        plt.show()
        




class Aminoacid():
    """
    
    """
    def __init__(self):
        """
        """
        self.location_x = None
        self.location_y = None
        self.step = None

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
    


protein = Protein("HPCHP")
protein.create_bonds()
# print(protein.sequence_list)

for acid in protein.sequence_list:
    print(f"Class = {acid}, Step = {acid.step}, Coordinates = ({acid.location_x}, {acid.location_y})")

##################################3



protein.visualize()
protein.create_output()

