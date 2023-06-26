from ..classes.protein import Protein

def UserInterface():
    print("Hello and welcome to our Protein folder!")
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

    algorithm = input("Which algorithm would you like to run?\n")
    
    #to do: Implement way to run algorithms

UserInterface()