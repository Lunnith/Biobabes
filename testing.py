from code.classes.protein import Protein
from code.visualisation.visualize import visualize

proteinA = Protein("HPCHPCHH")
proteinA.create_bonds()
# print(protein.sequence_list)

for acid in proteinA.sequence_list:
    print(f"Class = {acid}, Step = {acid.step}, Coordinates = ({acid.location_x}, {acid.location_y})")

proteinA.create_output()
visualize(proteinA)