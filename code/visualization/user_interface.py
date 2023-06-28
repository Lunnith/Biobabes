from ..classes.protein import Protein
from ..algorithms.randomise import random_reassignment
from ..algorithms.depth_first import DepthFirst
from ..algorithms.important_parts import ImportantParts
from ..algorithms.greedy import Greedy
from .visualize import *

class UserInterface():
    def __init__(self) -> None:
        print("Hello and welcome to our Protein folder!")
        self.run()
    
    def run(self):
        protein = self.determine_sequence()
        process_further = True

        while process_further:
            algorithm = self.determine_algorithm()

            if algorithm == 'a': 
                protein = self.use_random(protein)
            if algorithm == 'b': 
                protein = self.use_depth_first(protein)
            if algorithm == 'c': 
                protein = self.use_greedy(protein)
            if algorithm == 'd': 
                protein = self.use_hill_climber(protein)
            
            print(f"The lowest score that this algorithm found is {protein.score}.")
            asked_process = input("Do you wish to run another algorithm, continuing with this folded protein? \n\
                                  If yes, input 'Yes' or 'y'. Otherwise, the folded protein will be shown and the program will quit.").lower()
            if asked_process == "yes" or asked_process == 'y':
                continue
            else: process_further = False
        visualize_protein(protein, dimensions=3)    
        
    
    def determine_sequence(self):
        correct = False
        while not correct:
            own_sequence = input("Would you like to fold one of our given proteins, or would you like to fold your own sequence?\
                                For own sequence: Type your own sequence, for a given sequence, type 'N' or 'no'\n").upper()
            if own_sequence == "NO" or own_sequence == "N":
                sequence = input("Please, pick one of the following sequences:\n\
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
            print(f"The sequence to fold has a length of {len(sequence)}: {sequence}")
            correct = input("Is this correct? Please input 'Yes' or 'y'. ").lower()
            if correct == 'y' or correct == 'yes':
                correct = True
        return protein
    
    def determine_algorithm(self):
        continued = False

        while continued == False:
            print("The algorithms available are:\n\
                A) Random\n\
                B) Depth-First\n\
                C) Greedy\n\
                D) Hill-climber (with or without Simulated Annealing)")
            algorithm = input("For more information about an algorithm, or to run it, input either 'A', 'B', 'C' or 'D'.\n\"").lower()
            if algorithm == 'a': print("\nThis algorithm randomly assigns a direction for each bond of the protein. After a certain amount of iterations, it returns the protein with the best achieved score yet.")
            elif algorithm == 'b': print("\nThis algorithm computes every possible fold and is guaranteed to find the optimum folding. \nWARNING: This algorithm becomes VERY slow with a sequence length of more than 12.")
            elif algorithm == 'c': print("\nThis algorithm looks for the best possible next move, but greedy desicions at the beginning of the folding process might prevent it from getting an optimal score.")
            elif algorithm == 'd': print("\nThis algorithm makes small changes to the folding, hoping to improve the score. After a certain amount of iterations, it returns the best fold yet.")
            elif algorithm == 'q': break
            else: 
                print("That is not a valid option. If you wish to quit, input 'Q'.")

            verification = input("Do you wish to continue with this algorithm?\nIf yes, input 'Yes' or 'y'.").lower()
            if verification == 'yes' or verification == 'y':
                continued = True
            else: continued = False
        return algorithm
    
    def use_random(self, protein):
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
        return protein
    
    def use_depth_first(self, protein):
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
            return important_parts.protein
        elif pruning == 'P_pruning': depth_first.run(P_pruning=True)
        elif pruning == 'Directions_pruning': depth_first.run(directions_pruning=True)

        return depth_first.protein

    def use_greedy(self, protein):
        correct_kwargs = False

        while not correct_kwargs:
            splits_correct = False
            while not splits_correct:
                print("The Greedy algorithm by default uses a split size of 1. which means it only looks for the best score per addition of a single aminoacid. Increasing the split size increases the duration of the algorithm.")
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
        return greedy.protein

    def use_hill_climber(self, protein):
        return protein

    
UserInterface()