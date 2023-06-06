import matplotlib.pyplot as plt


def visualize(self):
    pos_x = []
    pos_y = []
    color = []
    for aminoacid in self.sequence_list:
        pos_x.append(aminoacid.location_x)
        pos_y.append(aminoacid.location_y)
        color.append(aminoacid.color)

    plt.plot(pos_x, pos_y, c = 'black', linestyle = '-')
    plt.scatter(pos_x, pos_y, c = color, markersize = 12)
    plt.grid(True, linestyle = '--')
    plt.show()