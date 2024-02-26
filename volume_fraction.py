#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 13:47:07 2024

@author: vova
"""
import numpy as np
from scipy.spatial import cKDTree


def gen_tree(filename):
    coordinates = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
                parts = line.strip().split()
                x = parts[3]
                y = parts[4]
                z = parts[5]
                coordinates.append([x,y,z])
                
    coordinates = np.array(coordinates, dtype=float)
        # Build a KD-tree from the freezed coordinates for efficient distance queries
    tree = cKDTree(coordinates)
    return tree



def volume_fraction(tree):

    x = np.random.uniform(-26, -9, 1000000)
    y = np.random.uniform(-6, 6, 1000000)
    z = np.random.uniform(-6, 6, 1000000)
    coords = np.column_stack([x,y,z]).reshape(-1, 3)
    
    distances, _ = tree.query(coords)
    
    mask = (distances >= 0.375)

    return np.mean(mask)

print(volume_fraction(gen_tree('linear_eq_coordinates.txt')))
print(volume_fraction(gen_tree('ring_eq_coordinates.txt')))
