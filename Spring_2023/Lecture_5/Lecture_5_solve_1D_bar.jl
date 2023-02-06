#EN.560.301 Structural Systems I 
#Lecture 5

using LinearAlgebra

#Calculate the axial displacement of a bar.

A = [1.0]
E = [29000.0]

node_coordinates = [[0.0, 0.0], [48.0, 0.0]]

element_connectivity = [(1, 2)]

supports = [1, 2, 4]

L = norm(node_coordinates[2] - node_coordinates[1])

