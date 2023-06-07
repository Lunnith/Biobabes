import matplotlib.pyplot as plt


def visualize(protein):
    """
    Visualize the folded protein with aminoacids and bonds.
    """
    pos_x = []
    pos_y = []
    color = []
    marker = []

    for aminoacid in protein.sequence_list:
        pos_x.append(aminoacid.location_x)
        pos_y.append(aminoacid.location_y)
        color.append(aminoacid.color)
        marker.append(aminoacid.type)

    plt.plot(pos_x, pos_y, c = 'black', linestyle = '-', linewidth = 0.7, zorder = 1)
    plt.scatter(pos_x, pos_y, c = color, m = marker, s = 50, zorder = 2)
    plt.grid(True, linestyle = '--')
    plt.show()