import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def visualize_protein2D(protein):
    """
    Visualize the folded protein with aminoacids and bonds.
    """
    pos_x = []
    pos_y = []

    for aminoacid in protein.sequence_list:
        pos_x.append(aminoacid.location[0])
        pos_y.append(aminoacid.location[1])
        plt.scatter(aminoacid.location[0], aminoacid.location[1], c = aminoacid.color, marker = f'${aminoacid.type}$', s = 200, zorder = 2)

    plt.plot(pos_x, pos_y, c = 'black', linestyle = '-', linewidth = 1.5, zorder = 1)
    visualize_bonds(protein)
    equal_axis(pos_x, pos_y)
    plt.grid(True, linestyle = '--')
    plt.show()

def equal_axis(x_list, y_list):
    """
    To_do: accept z_list. but also This function should be shortened.
    Make sure that the x-axis and y-axis are both of the same length, 
    so that the protein gets shown with equally long bonds.
    """
    max_x = max(x_list)
    max_y = max(y_list)
    min_x = min(x_list)
    min_y = min(y_list)
    length_x = max_x - min_x
    length_y = max_y - min_y
    
    if length_x < length_y:
        diff = length_y - length_x
        if diff % 2 == 0:
            max_x += diff // 2
            min_x -= diff // 2
        elif diff % 2 == 1:
            max_x += diff // 2 + 1
            min_x -= diff // 2

    elif length_y < length_x:
        diff = length_x - length_y
        if diff % 2 == 0:
            max_y += diff // 2
            min_y -= diff // 2
        elif diff % 2 == 1:
            max_y += diff // 2 + 1
            min_y -= diff // 2

    plt.xlim([min_x - 1, max_x + 1])
    plt.ylim([min_y - 1, max_y + 1])


def visualize_bonds(protein):
    """
    color hh = 'deepskyblue'
    color cc = 'lime'
    color ch = 'khaki'
    """
    for bond in protein.hh_ch_bonds:
        plt.plot((bond[0][0], bond[1][0]), (bond[0][1], bond[1][1]), linestyle='dotted', linewidth=2, c='slategrey')
    
    # for bond in protein.ch_bonds:
    #     plt.plot((bond[0][0], bond[1][0]), (bond[0][1], bond[1][1]), linestyle='dotted', linewidth=2, c='slategrey')

    for bond in protein.cc_bonds:
        plt.plot((bond[0][0], bond[1][0]), (bond[0][1], bond[1][1]), linestyle='dotted', linewidth=2, c='slategrey')


def visualize_scores(list_of_scores):
    sns.histplot(list_of_scores, kde=True)
    plt.xlabel("Score")
    plt.ylabel("Frequency")
    plt.show()

def visualize_protein(protein, dimensions):
    if dimensions == 2:
        visualize_protein2D(protein)
    elif dimensions == 3:
        visualize_protein3D(protein)
    else:
        print("ERROR: This function is written for either 2 or 3 dimensions.")
        return

def visualize_protein3D(protein):
    """
    To_do: implement bonds, improve visibility, see if you can make it interactive.
    Visualize the folded protein with aminoacids and bonds.
    """
    pos_x = []
    pos_y = []
    pos_z = []

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection = '3d')

    for aminoacid in protein.sequence_list:
        pos_x.append(aminoacid.location[0])
        pos_y.append(aminoacid.location[1])
        pos_z.append(aminoacid.location[2])
        ax.scatter(aminoacid.location[0], aminoacid.location[1], aminoacid.location[2], color = aminoacid.color)

    ax.plot(pos_x, pos_y, pos_z, color = 'black')
    fig.add_axes(ax)

    plt.show()