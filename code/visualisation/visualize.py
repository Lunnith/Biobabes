import matplotlib.pyplot as plt

def visualize_protein(protein):
    """
    Visualize the folded protein with aminoacids and bonds.
    """
    pos_x = []
    pos_y = []

    for aminoacid in protein.sequence_list:
        pos_x.append(aminoacid.location_x)
        pos_y.append(aminoacid.location_y)
        plt.scatter(aminoacid.location_x, aminoacid.location_y, c = aminoacid.color, marker = f'${aminoacid.type}$', s = 200, zorder = 2)

    plt.plot(pos_x, pos_y, c = 'black', linestyle = '-', linewidth = 1.5, zorder = 1)
    plt.grid(True, linestyle = '--')
    plt.show()

def visualise_bonds(protein):
    """
    
    """