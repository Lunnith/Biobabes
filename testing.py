from code.classes import protein
from code.classes import aminoacid
from code. visualisation import visualize

protein = Protein("HPCHP")
protein.create_bonds()
# print(protein.sequence_list)

for acid in protein.sequence_list:
    print(f"Class = {acid}, Step = {acid.step}, Coordinates = ({acid.location_x}, {acid.location_y})")

protein.visualize()
protein.create_output()