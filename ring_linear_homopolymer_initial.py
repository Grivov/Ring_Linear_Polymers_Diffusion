#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 12:06:11 2024
@author: vova
"""

import numpy as np
import matplotlib.pyplot as plt

def square_walk(length, diameter, x0, y0, z0):
    """
    Generate coordinates for a square walk.
    
    Parameters:
    length (int): The length of the walk.
    diameter (float): The diameter of the walk.
    x0, y0, z0 (float): Starting coordinates.

    Returns:
    tuple: Arrays of X, Y, Z coordinates.
    """
    half_side = length / 2 - 1
    x_frame = np.append(
        np.append(np.arange(0, half_side - 1, 1), [half_side - 1, half_side - 1]),
        np.append(np.arange(half_side - 1, 0, -1), [0, 0])
    )
    x = x_frame * diameter + x0

    y_frame = np.append(
        np.append(np.repeat(0, half_side - 1), [0, 1]),
        np.append(np.repeat(2, half_side - 1), [2, 1])
    )
    y = y_frame * diameter + y0

    z = np.repeat(z0, length)

    return x, y, z


def generate_coords(length, diameter, positions_y, positions_z):
    """
    Generate coordinates for multiple square walks.
    ...
    """
    half_side = length / 2 - 1
    x0_array = np.array([-diameter * half_side, 0.1]) + 0.3
    y0_array = diameter * positions_y + 0.4
    z0_array = diameter * positions_z + 0.45

    chain_number = 0
    coords_all = {'x': [], 'y': [], 'z': [], 'chain_number': []}

    for x0 in x0_array:
        for y0 in y0_array:
            for z0 in z0_array:
                chain_number += 1
                x, y, z = square_walk(length, diameter, x0, y0, z0)
                chain_number_list = [chain_number for _ in range(length)]

                coords_all['x'].append(x)
                coords_all['y'].append(y)
                coords_all['z'].append(z)
                coords_all['chain_number'].append(chain_number_list)

    for key in coords_all:
        coords_all[key] = np.array(coords_all[key]).flatten()

    plot_3d_scatter(coords_all)

    return coords_all


def plot_3d_scatter(coords):
    """
    Plot a 3D scatter plot of the coordinates.
    ...
    """
    fig = plt.figure(figsize=[10, 10])
    ax = fig.add_subplot(projection='3d')
    ax.scatter(
        coords['x'], coords['y'], coords['z'],
        s=10, c=coords['chain_number'], cmap='Set1', alpha=0.2
    )
    ax.plot3D(coords['x'], coords['y'], coords['z'], alpha=0.3)
    plt.show()


def create_file(name, mode, length, diameter, x_box_dim, positions_y, positions_z, linear):
    """
    Create a file with the generated coordinates.

    Parameters:
    name (str): The name of the file to create.
    mode (str): The file opening mode.
    length, diameter, x_box_dim (float): Parameters for coordinate generation.
    positions_y, positions_z (np.array): Arrays of positions.
    linear (bool): Flag to determine the structure type.
    """
    coords_all = generate_coords(length, diameter, positions_y, positions_z)

    polymers_num = 2 * len(positions_y) * len(positions_z)
    num_atoms = length * polymers_num
    num_bonds = num_atoms - polymers_num
    num_bonds_per_pol = length - 1

    box_dim_y = max(np.abs(coords_all['y'])) + 0.5 * diameter
    box_dim_z = max(np.abs(coords_all['z'])) + 0.5 * diameter

    with open(name, mode) as f:
        f.write(f'LAMMPS chain data file\n\n{num_atoms} atoms\n{num_bonds} bonds\n')
        f.write('0 angles\n0 dihedrals\n0 impropers\n\n1 atom types\n1 bond types\n')
        f.write('0 angle types\n0 dihedral types\n0 improper types\n\n')
        f.write(f'-{x_box_dim} {x_box_dim} xlo xhi\n')
        f.write(f'-{box_dim_y} {box_dim_y} ylo yhi\n')
        f.write(f'-{box_dim_z} {box_dim_z} zlo zhi\n\nMasses\n\n1 14.13717\n\nAtoms\n\n')

        for i in range(num_atoms):
            x, y, z = coords_all['x'][i], coords_all['y'][i], coords_all['z'][i]
            chain_number = coords_all['chain_number'][i]
            f.write(f'{i + 1} {chain_number} 1 {round(x, 4)} {round(y, 4)} {round(z, 4)}\n')

        f.write('\nVelocities\n\n')
        for i in range(num_atoms):
            f.write(f'{i + 1} 0 0 0\n')
        
        f.write('\nBonds\n\n')

                    
        if linear:
            bond_counter = 1
            for pol in range(polymers_num):
                for monomer in range(num_bonds_per_pol):
                        # Connect current monomer to the next monomer
                        first_atom = monomer + 1 + pol * num_bonds_per_pol + pol
                        second_atom = monomer + 2 + pol * num_bonds_per_pol + pol
                        f.write(f"{bond_counter} 1 {first_atom} {second_atom}\n")
                        bond_counter += 1
                    
        if not linear:
            num_bonds = num_atoms
            num_bonds_per_pol = length
            bond_counter = 1
            for pol in range(polymers_num):
                for monomer in range(num_bonds_per_pol):
                        # Connect current monomer to the next monomer
                        first_atom = monomer + 1 + pol * num_bonds_per_pol
                        second_atom = monomer + 2 + pol * num_bonds_per_pol if monomer < (num_bonds_per_pol - 1) else 1 + pol * num_bonds_per_pol
                        f.write(f"{bond_counter} 1 {first_atom} {second_atom}\n")
                        bond_counter += 1


# Main execution
name = 'test_initial_file.initial'
mode = 'w'
length = 36
diameter = 0.75
x_box_dim = 13 * 3
positions_y = np.arange(-12.2, 12.2, 3.05)
positions_z = np.arange(-12, 12, 1.2)
linear = False

create_file(name, mode, length, diameter, x_box_dim, positions_y, positions_z, linear)
