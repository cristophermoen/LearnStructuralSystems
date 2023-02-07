#EN.560.301 Structural Systems I 
#Lecture 5
#February 2023

using LinearAlgebra

#Calculate the axial displacement of a bar.

A = [1.0, 2.6, 4.5]
E = [29000.0, 4000.0]

node_coordinates = [[0.0, 0.0], [48.0, 0.0]]

element_connectivity = [(1, 2), (2,3), (3,4)]

supports = [1, 2, 4]

L = norm(node_coordinates[2] - node_coordinates[1])


