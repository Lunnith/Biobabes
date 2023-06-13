
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

def equal_axis(x_list, y_list, z_list=None):
    """
    To_do: This function should be shortened.
    Make sure that the x-axis, y-axis and z-axis are both of the same length, 
    so that the protein gets shown with equally long bonds.
    """
    max_x = max(x_list)
    max_y = max(y_list)
    min_x = min(x_list)
    min_y = min(y_list)
    length_x = max_x - min_x
    length_y = max_y - min_y
    max_length = max([length_x, length_y])

    if z_list != None:
        max_z = max(z_list)
        min_z = min(z_list)
        length_z = max_z - min_z
        max_length = max([max_length, length_z])
    

    if max_length == length_x:
        diff_xy = length_x - length_y
        if diff_xy % 2 == 0:
            max_y += diff_xy // 2
            min_y -= diff_xy // 2
        elif diff_xy % 2 == 1:
            max_y += diff_xy // 2 + 1
            min_y -= diff_xy // 2

        if z_list != None:
            diff_xz = length_x - length_z
            if diff_xz % 2 == 0:
                max_z += diff_xz // 2
                min_z -= diff_xz // 2
            elif diff_xz % 2 == 1:
                max_z += diff_xz // 2 + 1
                min_z -= diff_xz // 2 

    elif max_length == length_y:
        diff_xy = length_y - length_x
        if diff_xy % 2 == 0:
            max_x += diff_xy // 2
            min_x -= diff_xy // 2
        elif diff_xy % 2 == 1:
            max_x += diff_xy // 2 + 1
            min_x -= diff_xy // 2
        
        if z_list != None:
            diff_yz = length_y - length_z
            if diff_yz % 2 == 0:
                max_z += diff_yz // 2
                min_z -= diff_yz // 2
            elif diff_yz % 2 == 1:
                max_z += diff_yz // 2 + 1
                min_z -= diff_yz // 2  
    
    elif z_list != None and max_length == length_z:
        diff_xz = length_z - length_x
        if diff_xz % 2 == 0:
            max_x += diff_xz // 2
            min_x -= diff_xz // 2
        elif diff_xz % 2 == 1:
            max_x += diff_xz // 2 + 1
            min_x -= diff_xz // 2 

        diff_yz = length_z - length_y
        if diff_yz % 2 == 0:
            max_y += diff_yz // 2
            min_y -= diff_yz // 2
        elif diff_yz % 2 == 1:
            max_y += diff_yz // 2 + 1
            min_y -= diff_yz // 2  

    plt.xlim([min_x - 1, max_x + 1])
    plt.ylim([min_y - 1, max_y + 1])
    if z_list != None:
        return [min_z - 1, max_z + 1]



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
            axs.plot(xs, ys, zs, linestyle='dotted', linewidth=2, color='slategrey')
        for bond in protein.cc_bonds:
            xs = [bond[0][0], bond[1][0]]
            ys = [bond[0][1], bond[1][1]]
            zs = [bond[0][2], bond[1][2]]
            axs.plot(xs, ys, zs, linestyle='dotted', linewidth=2, color='slategrey')  
        return axs

    else:
        print("ERROR: This function is written for either 2 or 3 dimensions.")
        return


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
    ax.set_zlim3d(equal_axis(pos_x, pos_y, pos_z))
    fig.add_axes(ax).set_axis_off()
    fig.add_axes(visualize_bonds(protein, dimensions=3, axs=ax))
    plt.show()