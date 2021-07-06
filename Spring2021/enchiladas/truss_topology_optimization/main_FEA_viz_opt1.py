#!~/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 10:41:52 2018
Finite Element analysis in trusses with optimization
Author: @nikorose
"""
# =============================================================================
# Importing and Defining the libraries needed
# =============================================================================
import numpy as np
from FEAtrussDP3Opt import FEA #Vizualization
from scipy.optimize import minimize
import matplotlib.pyplot as plt


#Setting the classes needed
FEM = FEA()
#viz = Vizualization()

# =============================================================================
# Initial Variables statements
#  E: modulus of elasticity A: area of cross section L: length of bar 
# =============================================================================

# =============================================================================
# #Change the number of rows and columns for trusses
rows = 16
cols = 4
# =============================================================================
elementNodes, nodeCoordinates = FEM.structure2(rows,cols,1)
stress= np.zeros(elementNodes.shape[0])
# =============================================================================
#Change the Young Modulus here
E=1
# =============================================================================
# Defining areas for each element
A = np.ones((1,elementNodes.shape[0]))
# Defining element_length of each element
len_elements= np.array([[np.sqrt((e[0,0]-e[1,0])**2+(e[0,1]-e[1,1])**2) for e in elementNodes]] )
Vmax = np.dot(A,len_elements.T)

# =============================================================================
#  for structure:  displacements: displacement vector # force : force vector 
#  stiffness: stiffness matrix 
# =============================================================================

numberElements=elementNodes.shape[0]
numberNodes=nodeCoordinates.shape[0]
xx=nodeCoordinates[:,0]
yy=nodeCoordinates[:,1]

GDof=2*numberNodes
#U=np.zeros((GDof,1))
force=np.zeros((GDof,1))

''' ORIGINAL CODE
# Add applied load. Appoint gives the DOF it's applied to, here is 81)
appoint = (cols+1)*rows*2 + 1
force[appoint]=-1
# This defines the fixed DOF
prescribedDof=np.arange(cols*2+1)[np.newaxis].T;
index_elem =  FEM.indexation(elementNodes, nodeCoordinates)
'''

''' 
# BEAM PROBLEM
# Add applied load. DOF 49, negative direction
appoint = 49
force[appoint] = -10
# This defines the fixed DOF
prescribedDof = np.array([[0], [1], [81]])
index_elem =  FEM.indexation(elementNodes, nodeCoordinates)
'''

'''
#BEAM PROBLEM with 4x16
# Add applied load. DOF 89, negative direction
appoint = 89
force[appoint] = -3
# This defines the fixed DOF
prescribedDof = np.array([[0], [1], [161]])
index_elem =  FEM.indexation(elementNodes, nodeCoordinates)
'''

'''
#BEAM PROBLEM with 8x32
# Add applied load. DOF 255, negative direction
appoint = 255
force[appoint] = -1
# This defines the fixed DOF
prescribedDof = np.array([[0], [1], [511]])
index_elem =  FEM.indexation(elementNodes, nodeCoordinates)
'''

'''
#PLANE WING PROBLEM with 3x10
# Add applied load. DOF 89, negative direction
appoint = [7, 15, 23, 31, 39, 47, 55, 63, 71, 79 , 87]
force[appoint] = -.5
# This defines the fixed DOF
prescribedDof = np.array([[0], [1], [6], [7]])
index_elem =  FEM.indexation(elementNodes, nodeCoordinates)
'''

# BEAM PROBLEM with 4x16 DISTRIBUTED LOAD
# Add applied load. DOF 89, negative direction
appoint = 89
force[appoint] = -3
# This defines the fixed DOF
prescribedDof = np.array([[0], [1], [161]])
index_elem =  FEM.indexation(elementNodes, nodeCoordinates)

# =============================================================================
#  Optimization Section
# =============================================================================
def is_pos_def(x):
    '''
    To know if the matrix is positive definite.
    '''
    return np.all(np.linalg.eigvals(x) > 0)

def write(Xi):
    file = open('results.txt','a+') 
    file.write(str(Xi)+'\n') 
        
def optim_fun(x, *args):
    '''
    The Scipy Library is used to implement a solution for a 
    optimization problem looking for a compliant reduction.
    '''
    GDof, elementNodes, nodeCoordinates, indices, E, prescribedDof, force, writing = args
    stiff = FEM.formStiffness2Dtruss(x, GDof, elementNodes,nodeCoordinates, indices, E)
    displ = FEM.solution(GDof,prescribedDof,stiff, force)
    g0 = np.dot(force.T,displ)
    # If writing is activated, it will append the result of every iteration in the file results
    # you'll have to erase the file on every running  in order to not overwriting the last simulation
    if writing:  write(x)
    return g0[0][0]

# Defining the optimization problem and settings
cons = ({'type': 'ineq', 'fun': lambda x:  Vmax[0] - np.dot(x[:],len_elements.T)[0]}) 
# =============================================================================
# Change the bounds here
bnds = tuple([(.1,20)]*A.shape[1])
# =============================================================================

arguments = (GDof, elementNodes, nodeCoordinates, index_elem, E, prescribedDof, force, False)
options={'gtol': 1e-5, 'disp': True, 'maxiter':50 }


# =============================================================================
# Change the optimization method here
res = minimize(optim_fun, A, args=arguments, method='SLSQP', jac=None, 
               hess=None, hessp=None, bounds=bnds, constraints=cons, 
               tol=None, callback=None, options=options)
AreaOpt = res.x
Maxiter = res.nit
# =============================================================================

# =============================================================================
# boundary conditions and solution 
# =============================================================================

stiffness= FEM.formStiffness2Dtruss(AreaOpt,GDof,elementNodes, nodeCoordinates, \
                                    index_elem, E)

# solution
displacements = FEM.solution(GDof, prescribedDof, stiffness, force)

# stresses at elements
Stresses = FEM.stresses2Dtruss(elementNodes, index_elem, displacements, E)

# output displacements/reactions
Reactions = FEM.outputDisplacementsReactions(displacements,stiffness,GDof,prescribedDof)

# compliance Calculation
#Want to keep compliance between 1 and 100
g0 = np.dot(force.T,displacements)
print(g0)

# =============================================================================
# Plotting the behavior of the optimization result
# =============================================================================
#linewidth = 5.0
#scale = 2.0
#force_node = int((appoint-1)/2)
#prescribed_in_nodes = [int(e/2) for e in range(0,len(prescribedDof),2)]
#viz.visualizeAreas(elementNodes, nodeCoordinates, AreaOpt, prescribed_in_nodes,
                   #force_node, scale, Stresses, linewidth, Maxiter)

# =============================================================================
# Trying to display deformation and Areas
# =============================================================================


#viz.displayTruss(index_elem, nodeCoordinates, Stresses,Area=A, \
                 #name='Stress (MPa)', displace=True, displacements=displacements, \
                 #scale=0.1)

for i in range(0, len(AreaOpt)):
    x_vals = [elementNodes[i][0][0], elementNodes[i][1][0]]
    y_vals = [elementNodes[i][0][1], elementNodes[i][1][1]]

    color = 'blue'  # Blue for compression
    if Stresses[i][0] > 0.:
        color = 'red'

    plt.plot(x_vals, y_vals, color, linewidth = AreaOpt[i])

plt.show()

''' FOR DISPLAYING TWO PLOTS, STRESS AND AREA, SEPARATELY
figure, axis = plt.subplots(2, 1)

for i in range(0, len(AreaOpt)):
    x_vals = [elementNodes[i][0][0], elementNodes[i][1][0]]
    y_vals = [elementNodes[i][0][1], elementNodes[i][1][1]]
    axis[0].plot(x_vals, y_vals, 'gray', linewidth = AreaOpt[i])

axis[0].set_title("Optimized Areas")

for i in range(0, len(AreaOpt)):
    x_vals = [elementNodes[i][0][0], elementNodes[i][1][0]]
    y_vals = [elementNodes[i][0][1], elementNodes[i][1][1]]
    stress = Stresses[i][0]
    color = 'blue' #Blue for compression
    if stress > 0.:
        color = 'red' #Red for tension
    axis[1].plot(x_vals, y_vals, color, linewidth = 10 * stress)

axis[1].set_title("Stresses")

plt.show()
'''

print("hello")