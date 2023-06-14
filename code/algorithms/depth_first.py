import copy

class DepthFirst():
    """
    
    """
    def __init__(self, protein, dimensions):
        self.protein = protein

        self.directions = 
        self.states = [copy.deepcopy(self.protein)]

        self.best_state = None
        self.best_score = 0

    def get_next_state(self):
        """
        """
        return self.states.pop()
    
    def build_children(self, protein, aminoacid):
    
    def check_solution(self, new_fold):

    def run(self):
