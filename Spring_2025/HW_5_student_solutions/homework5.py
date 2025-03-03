# -*- coding: utf-8 -*-
"""
@author Emely Pacheco
Homework 5 - 3 bar truss analysis 
"""

import numpy as np
import matplotlib.pyplot as plt 

nodes = np.array([[0,0],[4,0],[2,3],])
elements = np.array([[1, 2, 1, 1],
                    [1, 3, 2, 2],
                    [2, 3, 3, 3]])
                   
material_properties = np.array([29000000.0, 10000.0, 1500000.0])
section_properties = np.array([1.0, 2.0, 3.0])


def truss_element(E, A, node1, node2):
    x1, y1 = nodes[node1 -1]
    x2, y2 = nodes[node2 - 1]
    
    L = np.sqrt((x2 -x1) ** 2 + (y2 - y1) ** 2)
    c = (x2 - x1) / L
    s = (y2 -y1) / L
    
    k_local = (E * A / L) * np.array([
        [c*c, c*s, -c*c, -c*s],
        [c*s, s*s, -c*s, -s*s],
        [-c*c, -c*s, c*c, c*s],
        [-c*s, -s*s, c*s, s*s]])
    
    return k_local, L, c, s


dof = len(nodes)* 2 
K_global = np.zeros((dof, dof))



for element_data in elements: 
    
    n1, n2, material_index, section_index = element_data.astype(int)
    E = material_properties[material_index -1]
    A = section_properties[section_index - 1]
    
    result = truss_element(E, A, n1, n2)
    k_local = result[0]
    
    dof2 = [ 2 * (n1 - 1), 2 * (n1 - 1) + 1,
           2 * (n2 - 1), 2 * (n2 - 1) + 1]
    
    for i in range(4):
        for j in range(4):
                K_global[dof2[i], dof2[j]] += k_local[i,j]
    
print("Global Stiffness Matrix of Truss")
print(K_global)
    

force_vector = np.zeros(dof)
fixed_dofs = [0,1]
other_dofs = []
for i in range(dof):
    if i not in fixed_dofs:
        other_dofs.append(i)

#attempting to solve for displacement
displacement = np.zeros(dof)
    
displacement[other_dofs] = np.linalg.solve(K_global[other_dofs][:, other_dofs], force_vector[other_dofs])


#attempting to solve for reaction forces 
reaction_f = []

for dof in fixed_dofs:
    reaction_force = 0
    
    for j in range(len(displacement)):
        reaction_f += K_global[dof, j] * displacement[j]

#attempting to solve for interal forces

for element_data in elements: 
    
    n1, n2, material_index, section_index = element_data.astype(int)
    E = material_properties[material_index -1]
    A = section_properties[section_index - 1]
    
    result = truss_element(E, A, n1, n2)
    k_local = result[0]
    
    dof2 = [ 2 * (n1 - 1), 2 * (n1 - 1) + 1,
           2 * (n2 - 1), 2 * (n2 - 1) + 1]
    L = result[1]
    
    internal_f = ((E*A) / L)
    
    

#plotting
plt.xlabel("X ")
plt.ylabel("Y ")
plt.title("Truss Element")
plt.plot(displacement[other_dofs], color='blue', linewidth = 2)
plt.show()




