from ..classes.protein import Protein
from ..algorithms.randomise import random_reassignment
from ..algorithms.depth_first import DepthFirst
from ..algorithms.important_parts import ImportantParts
from ..algorithms.greedy import Greedy
from ..algorithms.hill_climber import Hill_climber
from ..algorithms.simulated_annealing import SimulatedAnnealing
from .visualize import *
import copy


class UserInterface():
    """
    This class is defined to automatically run the userinterface for the protein folding program.
    """
    def __init__(self) -> None:
        print("Hello and welcome to our Protein folder! (At any time if you want to quit, just insert 'Q'!)")
        self.folded = False
        self.run()
    
    def run(self):
        """
        Run the user interface.
        """
        original_protein = self.determine_sequence()
        protein = copy.deepcopy(original_protein)
        process_further = True

        while process_further:
            if not self.folded: algorithm = self.determine_algorithm_unfolded()
            elif self.folded: algorithm = self.determine_algorithm_folded(protein)

            if algorithm == 'refolded':
                protein = copy.deepcopy(original_protein)
                continue
            if algorithm == 'stop': break
            protein = self.run_algorithm(algorithm, protein)
            
            print(f"\nThe lowest score that this algorithm found is {protein.score}.")
            asked_process = input("Do you wish to run another algorithm, continuing with this folded protein? \n\
                                  If yes, input 'Yes' or 'y'. Otherwise, the folded protein will be shown and the program will quit.\n").lower()
            if asked_process == "yes" or asked_process == 'y':
                continue
            else: process_further = False

        if self.folded: visualize_protein(protein, dimensions=3)
        else: print("\nThis protein is not folded and thus cannot be visualised.")  

    def run_algorithm(self, algorithm, protein): 
        """
        Initiate the start of the chosen algorithm
        """
        if algorithm == 'a': 
            protein = self.use_random(protein)
        if algorithm == 'b': 
            protein = self.use_depth_first(protein)
        if algorithm == 'c': 
            protein = self.use_greedy(protein)
        if algorithm == 'd': 
            protein = self.use_hill_climber(protein)   
        return protein  
    
    def determine_sequence(self):
        """
        Ask the user what sequence to use.
        """
        correct = False
        while not correct:
            print("     Would you like to fold one of our given proteins? Insert 'Yes' or 'y'. ")
            own_sequence = input("     If you would like to fold your own sequence, please insert your sequence: \n").upper()
            if own_sequence == "YES" or own_sequence == "Y":
                sequence = input("\nPlease, pick one of the following sequences:\n\
                Sequences excluding C's: \n\
                A) HHPHHHPH \n\
                B) HHPHHHPHPHHHPH \n\
                C) HPHPPHHPHPPHPHHPPHPH \n\
                D) PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP \n\
                E) HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH \n\
                Sequences including C's: \n\
                F) PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP \n\
                G) CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC \n\
                H) HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH \n\
                I) HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH \n").lower()
                if sequence == 'a': sequence = 'HHPHHHPH'
                elif sequence == 'b': sequence = 'HHPHHHPHPHHHPH'
                elif sequence == 'c': sequence = 'HPHPPHHPHPPHPHHPPHPH'
                elif sequence == 'd': sequence = 'PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP'
                elif sequence == 'e': sequence = 'HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH'
                elif sequence == 'f': sequence = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'
                elif sequence == 'g': sequence = 'CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC'
                elif sequence == 'h': sequence = 'HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH'
                elif sequence == 'i': sequence = 'HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH'
                
            else:
                sequence = own_sequence
            
            protein = Protein(sequence)
            print(f"\nThe sequence to fold has a length of {len(sequence)}: {sequence}")
            correct = input("Is this correct? Please input 'Yes' or 'y'. \n").lower()
            if correct == 'y' or correct == 'yes':
                correct = True
        return protein
    

    def determine_algorithm_unfolded(self):
        """
        Only gets called upon if self.folded = False.
        Asks the user which algorithm to run.
        """
        continued = False

        while continued == False:
            print("\nThe algorithms available are:\n\
                A) Random\n\
                B) Depth-First\n\
                C) Greedy\n\
                D) Hill-climber (with or without Simulated Annealing)")
            algorithm = input("For more information about an algorithm, or to run it, input either 'A', 'B', 'C' or 'D'.\n\"").lower()
            if algorithm == 'q': break
            elif algorithm == 'a' or algorithm == 'b' or algorithm == 'c' or algorithm == 'd': continued=True
            else: continued = False
        return algorithm
    
    def determine_algorithm_folded(self, protein):
        """
        Only gets called upon if self.folded = True
        Asks the user how to continue.
        """
        print("\n\nAs the protein is already folded, not every algorithm can be run anymore.")
        print("The only algorithm that can improve the pre-folded protein is the Hill-Climber (with, or without Simulated Annealing) algorithm.")
        print("However, if the score you have found is too low and you still want to try the other algorithms, we can unfold the protein for you.")
        choice_correct = False
        while not choice_correct:
            choice = input("\nDo you wish to unfold your protein and run another algorithm? If yes, input 'Yes' or 'y'.\n\
                        If you want to run the Hill-Climber algorithm, input 'Run' or 'r'. \n\
                        If you don't want to run the Hill-Climber and keep the folding like this, input 'No' or 'n'.\n").lower()
            if choice == 'yes' or choice == 'y':
                choice_correct = True
                choice_sure = input(f"Are you sure that you want to reject the fold (score={protein.score}) you have gotten until now? ").lower()
                if choice_sure == 'yes' or choice_sure == 'y':
                    self.folded = False
                    algorithm = 'refolded'
                    return algorithm
            elif choice == "run" or choice == 'r': algorithm = 'd'
            elif choice == "no" or choice == 'n' or choice == 'q': algorithm = 'stop'

            return algorithm


    def use_random(self, protein):
        """
        Asks the user for the parameters nessecary and runs the Random algorithm.
        """
        print("\nThis algorithm randomly assigns a direction for each bond of the protein. After a certain amount of iterations, it returns the protein with the best achieved score yet.")
        correct_kwargs = False
        while not correct_kwargs:
            print("\nThe only parameter to give the random algorithm is the amount of iterations you would like to run.")
            kwargs = input("How many iterations would you like to use? \n").lower()
            if kwargs == 'q': return protein
            try: kwargs = int(kwargs)
            except: continue
            else:
                verification = input(f"Now using {kwargs} iterations, do you wish to continue?\n").lower()
                if verification == "yes" or verification == 'y':
                    correct_kwargs = True
        protein = random_reassignment(protein, 3, k=kwargs)[0]
        self.folded = True
        return protein
    
    def use_depth_first(self, protein):
        """
        Asks the user for the parameters nessecary and runs the Depth-First algorithm.
        """
        print("\nThis algorithm computes every possible fold and is guaranteed to find the optimum folding. \nWARNING: This algorithm becomes VERY slow with a sequence length of more than 12.")
        depth_first = DepthFirst(protein, 3)
        important_parts = ImportantParts(protein, 3)

        correct_kwargs = False
        while not correct_kwargs:
            pruning_correct = False
            while not pruning_correct:
                print("\nThe depth first algorithm has serveral options:")
                print("For example, the sequence can automatically be split into important parts, to be able to handle longer sequences.")
                important_parts_answer = input("Would you like to use this option? If yes, type 'Yes' or 'y'.\n").lower()
                if important_parts_answer == 'yes' or important_parts_answer == 'y':
                    pruning = 'Important_parts'
                    pruning_correct = True
                    continue

                print("This algorithm can use P_pruning or Directions_pruning, which makes it only a little bit faster.")
                kwargs = input("Do you want to use pruning?\n").lower()
                if kwargs == 'q': return protein

                elif kwargs == 'yes' or kwargs == 'y': 
                    pruning = input("Do you want to use A) P_pruning, or B) Directions_pruning?\n").lower()
                    if pruning == 'a' or pruning == 'b': pruning_correct = True
                    continue
                
                elif kwargs == 'no' or 'n': 
                    pruning = None
                    pruning_correct = True
                    continue

            if pruning == 'q': return protein
            elif pruning == 'a': pruning = 'P_pruning'
            elif pruning == 'b': pruning = 'Directions_pruning'
            
            kwargs = input(f"Now using the following pruning method: {pruning}. Correct?\n").lower()
            if kwargs == "yes" or kwargs == "y": correct_kwargs = True

        if pruning == None: depth_first.run()
        elif pruning == 'Important_parts': 
            important_parts.run()
            self.folded = True
            return important_parts.protein
        elif pruning == 'P_pruning': depth_first.run(P_pruning=True)
        elif pruning == 'Directions_pruning': depth_first.run(directions_pruning=True)

        self.folded = True
        return depth_first.protein

    def use_greedy(self, protein):
        """
        Asks the user for the parameters nessecary and runs the Greedy algorithm.
        """
        print("\nThis algorithm looks for the best possible next move, but greedy desicions at the beginning of the folding process might prevent it from getting an optimal score.")
        correct_kwargs = False

        while not correct_kwargs:
            splits_correct = False
            while not splits_correct:
                print("\nThe Greedy algorithm by default uses a split size of 1. which means it only looks for the best score per addition of a single aminoacid. Increasing the split size increases the duration of the algorithm.")
                splits = input("Do you want to increase the split size? If yes, input 'Yes' or 'y'.\n").lower()
                if splits == 'q': return protein
                if splits == "yes" or splits == "y":
                    split_size = int(input("Please input the desired split size: "))
                    splits_correct = True
                    break
                else: 
                    split_size = 1
                    splits_correct=True
                    break
            
            kwargs = input(f"Now using a split size of {split_size}, correct?\n").lower()
            if kwargs == "yes" or kwargs == "y":
                splits_correct = True
                correct_kwargs = True

        greedy = Greedy(protein, 3, splits=split_size)
        greedy.run()
        self.folded = True
        return greedy.protein

    def use_hill_climber(self, protein):
        """
        Asks the user for the parameters nessecary and runs the Hill Climber algorithm.
        """
        print("\nThis algorithm makes small changes to the folding, hoping to improve the score. After a certain amount of iterations, it returns the best fold yet.")
        simanneal = False
        correct_kwargs = False
        folded = self.folded

        hill_climber = Hill_climber(protein, folded=folded)
        print("The Hill-Climber algorithm without Simulated Annealing rejects all changes that have a negative impact. Because of this, the Hill-Climber could end up at a local optimum. The version with Simulated Annealing sometimes accepts a worse score, in order to be able to get out of these local optima and hopefully end up in the global optimum.")
        
        while not correct_kwargs:
            anneal = input("would you like to use the Simulated Annealing addition? If yes, input 'Yes' or 'y'.\n").lower()
            if anneal == 'q': return protein
            if anneal == "yes" or anneal == 'y':
                simanneal = True
            
            try: n = int(input("How many bonds would you like to change at the same time?\n"))
            except:
                print("Please insert a number.")
                continue
            if n > 10:
                print("Please insert a maximum of 10")
                continue

            try: iterations = int(input("How many iterations do you want to execute?\n"))
            except:
                print("Please insert a number.")
                continue

            kwargs = input(f"Now using {iterations} iterations, changing {n} bonds, and the use of Simulated Annealing is {simanneal}, correct?\n").lower()
            if kwargs == "yes" or kwargs == "y":
                correct_kwargs = True
            else: simanneal = False
        
        if simanneal:
            hill_climber = SimulatedAnnealing(protein, start_n=n, folded=folded)
        protein = hill_climber.run_i_iterations(protein, iterations, n, sim_annealing=simanneal)[0]
        self.folded = True
        return protein