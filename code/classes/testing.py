from . import aminoacid
from . import protein

protein = Protein("HPCHP")
protein.create_bonds()
# print(protein.sequence_list)

for acid in protein.sequence_list:
    print(f"Class = {acid}, Step = {acid.step}, Coordinates = ({acid.location_x}, {acid.location_y})")

##################################3



protein.visualize()
protein.create_output()

