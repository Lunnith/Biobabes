
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def visualize_protein(protein, dimensions):
    if dimensions == 2:
        visualize_protein2D(protein)
    elif dimensions == 3:
        visualize_protein3D(protein)
    else:
        print("ERROR: This function is written for either 2 or 3 dimensions.")
        return

def visualize_protein2D(protein):
    """
    Visualize the folded protein with aminoacids and bonds in 2D.
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


def visualize_protein3D(protein):
    """
    Visualize the folded protein with aminoacids and bonds in 3D.
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
        ax.scatter(aminoacid.location[0], aminoacid.location[1], aminoacid.location[2], color = aminoacid.color, marker = f'${aminoacid.type}$', s=50)

    ax.plot(pos_x, pos_y, pos_z, color = 'black', linewidth=1.5)
    ax.set_zlim3d(equal_axis(pos_x, pos_y, pos_z))
    fig.add_axes(ax).set_axis_off()
    fig.add_axes(visualize_bonds(protein, dimensions=3, axs=ax))
    plt.show()


def equal_axis(x_list, y_list, z_list=None):
    """
    Make sure that the x-axis, y-axis and z-axis are both of the same length, 
    so that the protein gets shown with equally long bonds.
    """
    min_x, max_x, length_x = get_min_max_length(x_list)
    min_y, max_y, length_y = get_min_max_length(y_list)
    max_length = max([length_x, length_y])
    if z_list != None:
        min_z, max_z, length_z = get_min_max_length(z_list)
        max_length = max([max_length, length_z])

    if max_length == length_x:
        diff_xy = length_x - length_y
        min_y, max_y = calculate_axis_change(diff_xy, min_y, max_y)

        if z_list != None:
            diff_xz = length_x - length_z
            min_z, max_z = calculate_axis_change(diff_xz, min_z, max_z)

    elif max_length == length_y:
        diff_xy = length_y - length_x
        min_x, max_x = calculate_axis_change(diff_xy, min_x, max_x)
        
        if z_list != None:
            diff_yz = length_y - length_z
            min_z, max_z = calculate_axis_change(diff_yz, min_z, max_z)
    
    elif z_list != None and max_length == length_z:
        diff_xz = length_z - length_x
        min_x, max_x = calculate_axis_change(diff_xz, min_x, max_x) 

        diff_yz = length_z - length_y
        min_y, max_y = calculate_axis_change(diff_yz, min_y, max_y) 

    plt.xlim([min_x - 1, max_x + 1])
    plt.ylim([min_y - 1, max_y + 1])
    if z_list != None:
        return [min_z - 1, max_z + 1]
        
def get_min_max_length(list):
    min_list = min(list)
    max_list = max(list)
    return min_list, max_list, (max_list - min_list)

def calculate_axis_change(diff, min_smaller_axis, max_smaller_axis):
        if diff % 2 == 0:
            max_smaller_axis += diff // 2
            min_smaller_axis -= diff // 2
        elif diff % 2 == 1:
            max_smaller_axis += diff // 2 + 1
            min_smaller_axis -= diff // 2
        return min_smaller_axis, max_smaller_axis


def visualize_bonds(protein, dimensions=2, axs=None):
    """
    color hh = 'deepskyblue'
    color cc = 'lime'
    color ch = 'khaki'
    """
    if dimensions == 2:
        for bond in protein.hh_ch_bonds:
            plt.plot((bond[0][0], bond[1][0]), (bond[0][1], bond[1][1]), linestyle='dotted', linewidth=2, c='slategrey')

        for bond in protein.cc_bonds:
            plt.plot((bond[0][0], bond[1][0]), (bond[0][1], bond[1][1]), linestyle='dotted', linewidth=2, c='slategrey')

    elif dimensions == 3:
        for bond in protein.hh_ch_bonds:
            xs = [bond[0][0], bond[1][0]]
            ys = [bond[0][1], bond[1][1]]
            zs = [bond[0][2], bond[1][2]]
            axs.plot(xs, ys, zs, linestyle='dotted', linewidth=1, color='slategrey')
        for bond in protein.cc_bonds:
            xs = [bond[0][0], bond[1][0]]
            ys = [bond[0][1], bond[1][1]]
            zs = [bond[0][2], bond[1][2]]
            axs.plot(xs, ys, zs, linestyle='dotted', linewidth=1, color='slategrey')  
        return axs

    else:
        print("ERROR: This function is written for either 2 or 3 dimensions.")
        return


def visualize_scores(list_of_scores):
    sns.histplot(list_of_scores, kde=True)
    plt.xlabel("Score")
    plt.ylabel("Frequency")
    plt.show()