from protein import Protein


protein = Protein("HPCHP")
protein.create_bonds()
# print(protein.sequence_list)
protein.check_interactions(protein.sequence_list[3])

for acid in protein.sequence_list:
    print(f"Class = {acid}, Step = {acid.step}, Coordinates = ({acid.location_x}, {acid.location_y})")

protein.visualize()
protein.create_output()



