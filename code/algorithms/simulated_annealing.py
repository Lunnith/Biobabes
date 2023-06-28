from .hill_climber import Hill_climber
from ..classes.protein import Protein

import random
import math

class SimulatedAnnealing(Hill_climber):
    """
    A SimulatedAnnealing class to change the direction of one or bonds at a time to 
    a random other direction. Each improvement or equivalent solution is kept for the
    next iteration, sometimes it accepts worse solutions, depending on the current temperature.
    
    Most of the functions are similar to those of the HillClimber class, which is why
    we use that as a parent class.
    ...
    Attributes:
    -----------
    Inherits all attributes from Hillclimber. Extra attributes:
    T0: int
        starting temperature
    T: int
        current temperature
    temp_scheme: str
        cooling scheme
    start_n: int
        amount of bonds to change to start with

    Methods:
    -----------
    Inherits all methods from DepthFirst. Extra methods:
    update_temperature(alpha):
        implements cooling scheme and decreases number of bonds to be changed
    check_solution(new_folding):
        checks and accepts new foldings based on temperature
    reset_temperature(new_t, new_n, reset_protein):
        resets temperature, usefull during experiments
    """
    def __init__(self, protein: Protein, start_n: int, folded: bool = False, dimensions: int = 3, temperature: int = 1, prints: bool = False, temp_scheme: str = 'exp') -> None:
        # use init of Hillclimber class
        super().__init__(protein, dimensions, folded=folded, prints=prints)

        # starting temperature and current temperature
        self.T0 = temperature
        self.T = temperature
        self.temp_scheme = temp_scheme

        self.start_n = start_n
    
    def update_temperature(self, alpha=0.99) -> None:
        """
        Method to implement a linear or exponential cooling scheme. Temperature becomes zero after all iterations
        passed to the run() method have passed. Also smoothly decreases the number of bonds to be changed per iterations while temperature decreases.
        """
        # decrease the number of changed bonds when the temperature decreases
        if self.T0 == self.T:
            self.n_float = self.start_n

        if self.n > 1:
            beta = self.start_n / self.iterations

            self.n_float = self.n_float - beta
            self.n = round(self.n_float)
            
        # implements exponential or linear scheme
        if self.temp_scheme == 'exp':
            if self.T > 0.01:
                self.T = self.T * alpha
            else:
                self.T = 0.01

        elif self.temp_scheme == 'lin':
            alpha = self.T0 / self.iterations

            if self.T > 0.01:
                self.T = self.T - alpha
            else:
                self.T = 0.01

    def check_solution(self, new_folding: Protein) -> None:
        """
        Checks and accepts better solutions than the current solution.
        Sometimes accepts worse solutions, which depends on the current temperature.
        """
        new_score = new_folding.score
        old_score = self.protein.score

        # calculate probability of accepting new folding
        delta = -old_score + new_score
        try:
            probability = math.exp(-delta / self.T)
        except:
            probability = 1
 
        # pull a random number between 0 and 1 and see if we accept the graph
        if random.random() < probability:
            self.protein = new_folding
            self.lowest_score = new_folding.score
            self.improvement.append("Y") #Change got accepted
            
            if self.prints: print(f"Score updated to", new_score, end=" ")
            updated = True
        else: 
            self.improvement.append("N")
            updated = False
        
        # update temperature
        self.update_temperature()
        return updated

    def reset_temperature(self, new_t=None, new_n=None, reset_protein=None):
        """
        This can be run when experimenting, it resets the temperature after a run back to the start_temperature.
        It has also got a function to manipulate the temperature or the start_n for the next run.
        """
        if new_t is not None:
            self.T = new_t
        else: self.T = self.T0

        if new_n is not None:
            self.start_n = new_n
        
        if reset_protein is not None:
            self.protein = reset_protein