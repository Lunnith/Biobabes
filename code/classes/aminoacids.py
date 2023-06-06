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
        
        for acid in self.sequence_list:
            acid.step = random.choice(self.directions)

            if acid == self.sequence_list[0]:
                pass

            else: 
                step_taken = previous_acid.step
                old_location_x = previous_acid.location_x
                old_location_y = previous_acid.location_y

                if step_taken == -1:
                    acid.location_x = old_location_x - 1
                    acid.location_y = old_location_y
            
                elif step_taken == 1:
                    acid.location_x = old_location_x + 1
                    acid.location_y = old_location_y

                elif step_taken == -2:
                    acid.location_x = old_location_x
                    acid.location_y = old_location_y - 1

                elif step_taken == 2:
                    acid.location_x = old_location_x
                    acid.location_y = old_location_y + 1

            previous_acid = acid



    def check_bonds(self):
        """
        """
        pass

    def create_output(self):
        """
        """
    
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
# print(protein.sequence_list)

for acid in protein.sequence_list:
    print(f"Class = {acid}, Step = {acid.step}, Coordinates = ({acid.location_x}, {acid.location_y})")

##################################3



protein.visualize()


