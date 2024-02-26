#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:52:35 2024

@author: vova
"""
import numpy as np

def gen_coords(filename):
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
    return coordinates



def number_density(coordinates):

    condition_x = (coordinates[:, 0] > -25) & (coordinates[:, 0] < -9)
    condition_y = (coordinates[:, 1] > -5) & (coordinates[:, 1] < 5)
    condition_z = (coordinates[:, 2] > -5) & (coordinates[:, 2] < 5)
    
    # Combine conditions using logical AND
    mask = condition_x & condition_y & condition_z

    return sum(mask) / 16 / 10 / 10

print(number_density(gen_coords('/Users/vova/Downloads/linear_eq_coordinates.txt')))
print(number_density(gen_coords('/Users/vova/Downloads/ring_eq_coordinates.txt')))
