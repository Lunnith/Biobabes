import random
import math

from .hill_climber import Hill_climber
from ..classes.protein import Protein

class SimulatedAnnealing(Hill_climber):
    """
    A SimulatedAnnealing class to change the direction of one or bonds at a time to 
    a random other direction. Each improvement or equivalent solution is kept for the
    next iteration, sometimes it accepts worse solutions, depending on the current temperature.
    
    Most of the functions are similar to those of the HillClimber class, which is why
    we use that as a parent class.
    """
    def __init__(self, protein: Protein, start_n: int, folded: bool = False, dimensions: int = 3, temperature: int = 1) -> None:
        # use init of Hillclimber class
        super().__init__(protein, dimensions, folded=folded)

        # starting temperature and current temperature
        self.T0 = temperature
        self.T = temperature

        self.start_n = start_n
    
    def update_temperature(self) -> None:
        """
        Method to implement a linear cooling scheme. Temperature becomes zero after all iterations
        passed to the run() method have passed.
        """
        if self.T0 == self.T:
            self.n_float = self.start_n

        if self.n > 1:
            beta = self.start_n / self.iterations

            self.n_float = self.n_float - beta
            self.n = round(self.n_float)
            
        alpha = 0.99
        if self.T > 0.01:
            self.T = self.T * alpha
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
        probability = math.exp(-delta / self.T)
 
        # pull a random number between 0 and 1 and see if we accept the graph
        if random.random() < probability:
            self.protein = new_folding
        
        # update temperature
        self.update_temperature()
