import matplotlib.pyplot as plt


def visualize(self):
    pos_x = []
    pos_y = []
    color = []
    for aminoacid in self.sequence:
        pos_x.append(aminoacid.pos_x)
        pos_y.append(aminoacid.pos_y)
        color.append(aminoacid.color)

    plt.scatter(pos_x, pos_y, color)
    plt.plot(pos_x, pos_y, 'black-')
    plt.grid(True)
    plt.show()