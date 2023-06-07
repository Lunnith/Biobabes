from Aminoacid import *

import matplotlib.pyplot as plt
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
        self.sequence_list = []
        self.score = 0

        self.add_aminoacid(sequence)
        self.initialize_neighbours()


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
    
    def initialize_neighbours(self):
        """
        
        """
        for i in range(len(self.sequence_list)):
            self.sequence_list[i].neighbour1 = self.sequence_list[i - 1]
            if i < len(self.sequence_list) - 1:
                self.sequence_list[i].neighbour2 = self.sequence_list[i + 1]

        self.sequence_list[0].neighbour1 = None

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
                previous_acid.step = direction[2]

                # Set location of the acid
                acid.location_x = previous_acid.location_x + direction_x
                acid.location_y = previous_acid.location_y + direction_y


            previous_acid = acid

   
    def check_interactions(self, aminoacid):
        """
        Check surrounding of aminoacid for other aminoacids, if they are present, check type
        and change score according to the interaction type.
        """

        for acid in self.sequence_list:
            if acid != aminoacid.neighbour1 or acid != aminoacid.neighbour2:

                if acid.distance(aminoacid) == 1:
                    potential_interactor = acid

                if type(aminoacid) == Hydrophobic() and type(potential_interactor) == Hydrophobic():
                    self.score -= 1
        
                if type(aminoacid) == Hydrophobic() and type(potential_interactor) == Cysteine():
                    self.score -= 1
        
                if type(aminoacid) == Cysteine() and type(potential_interactor) == Cysteine():
                    self.score -= 5
            
            potential_interactor = None

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

        plt.plot(pos_x, pos_y, c = 'black', linestyle = '-', linewidth = 0.7, zorder = 1)
        plt.scatter(pos_x, pos_y, c = color, s = 50, zorder = 2)
        plt.grid(True, linestyle = '--')
        plt.show()
        