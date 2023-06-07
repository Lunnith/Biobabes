from code.classes import protein
from code.visualisation import visualize

proteinA = Protein("HPCHP")
proteinA.create_bonds()
# print(protein.sequence_list)

for acid in proteinA.sequence_list:
    print(f"Class = {acid}, Step = {acid.step}, Coordinates = ({acid.location_x}, {acid.location_y})")

proteinA.create_output()
visualize(proteinA)