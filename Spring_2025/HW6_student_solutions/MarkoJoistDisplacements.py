#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 13:02:24 2025

@author: evankelly
"""

import math
import numpy as np
import pandas as pd

#define truss nodes
nodesraw = pd.read_csv("nodes.csv")

elementsraw = pd.read_csv("elements.csv")

nodes = nodesraw.values

elements = elementsraw.values


#E in PSI
E = 29000000.0

#cross-sectional area A in in^2
A = [0.877, 0.30]

#calculate element length
l = np.zeros(43)
for i in range(len(elements)):
    
   x1 = nodes[elements[i][1]-1][1]
   x2 = nodes[elements[i][2]-1][1]
   y1 = nodes[elements[i][1]-1][2]
   y2 = nodes[elements[i][2]-1][2]
   
   length = math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2))
   
   l[i] = length
    

#calculate element orientation
theta = np.zeros(43)
for i in range(len(elements)):
    
   x1 = nodes[elements[i][1]-1][1]
   x2 = nodes[elements[i][2]-1][1]
   y1 = nodes[elements[i][1]-1][2]
   y2 = nodes[elements[i][2]-1][2]
   
   delta = [x2 - x1, y2 - y1]
   
   t = math.atan2(delta[1], delta[0])
   
   theta[i] = t
   

#calculate ke local
ke_local = []
for i in range(len(elements)):
    
    ke_local.append(E * A[elements[i][4]-1] * l[i])
    ke_local.append(0)
    ke_local.append(-E * A[elements[i][4]-1] * l[i])
    
    for j in range(5):
        ke_local.append(0)
        
    ke_local.append(-E * A[elements[i][4]-1] * l[i])
    ke_local.append(0)
    ke_local.append(E * A[elements[i][4]-1] * l[i])
    
    for k in range(5):
        ke_local.append(0)
    
ke_local = np.array(ke_local)

ke_local = ke_local.reshape((43,4,4))

#define element rotation matrix
beta = []
for i in range(len(elements)):
    
    beta.append(math.cos(theta[i]))
    beta.append(math.sin(theta[i]))
    
    for j in range(2):
        beta.append(0)
        
    beta.append(-math.sin(theta[i]))
    beta.append(math.cos(theta[i]))
    
    for k in range(4):
        beta.append(0)

    beta.append(math.cos(theta[i]))
    beta.append(math.sin(theta[i]))
    
    for m in range(2):
        beta.append(0)

    beta.append(-math.sin(theta[i]))
    beta.append(math.cos(theta[i]))
    
beta = np.array(beta)

beta = beta.reshape((43,4,4))


#set up degrees of freedom and big mac
dof_i = []
dof_j = []
for i in range(1, len(nodes)+1):
    dof_i.append(2*i - 1)
    dof_j.append(2*i)

dof = sorted(dof_i + dof_j)

num_dof = len(dof)


elements_dof = []
for i in range(len(elements)):
    
    elements_dof.append(dof_i[elements[i][1] - 1])
    elements_dof.append(dof_j[elements[i][1] - 1])
    elements_dof.append(dof_i[elements[i][2] - 1])
    elements_dof.append(dof_j[elements[i][2] - 1])
    
elements_dof = np.array(elements_dof)

elements_dof = elements_dof.reshape(43,4)

K = np.zeros((num_dof, num_dof))

#define global element stiffness matrix and combine into big mac
for i in range(len(elements)):
    
    betaT = np.transpose(beta[i])
    
    ke_global = betaT @ ke_local[i] @ beta[i]
    
    for v1 in range(4):
        for v2 in range(4):
            
            K[elements_dof[i][v1] - 1, elements_dof[i][v2] - 1] += ke_global[v1, v2]
            

#partition big mac and define external forces

w = 300/12

F = np.zeros(46)

top_chord = list(range(0,12))

for i in top_chord:
    
    F[2*i-1] += w * l[i-1]/2
    
    F[2*(i+1)-1] += w * l[i-1]/2

fixedDOF = [0, 1, 23]

freeDOF = []

for i in range(46):
    freeDOF.append(i)

for num in fixedDOF:
    freeDOF.remove(num)
    
Kpp1 = K[2:23, 2:23]

Kpp2 = K[24:46, 2:23]

Kpp3 = K[2:23, 24:46]

Kpp4 = K[24:46, 24:46]

Kpp12 = np.concatenate((Kpp1, Kpp2))

Kpp34 = np.concatenate((Kpp3, Kpp4))

Kpp = np.concatenate((Kpp12, Kpp34), axis=1)



Kss1 = K[0:2, 0:2]

Kss2 = K[23, 0:2]

Kss3 = K[0:2, 23]

Kss4 = K[23, 23]

Kss2 = Kss2.reshape((2,1))

Kss12 = np.append(Kss1, Kss2, axis=1)

Kss34 = np.append(Kss3, Kss4)

Kss34 = Kss34.reshape((1,3))

Kss = np.concatenate((Kss12, Kss34), axis=0)


Ksp1 = K[0:2, 2:23]

Ksp2 = K[0:2, 24:46]

Ksp3 = K[23, 2:23]

Ksp4 = K[23, 24:46]

Ksp3 = Ksp3.reshape((1,21))

Ksp4 = Ksp4.reshape((1,22))

Ksp12 = np.append(Ksp1, Ksp2, axis=1)

Ksp34 = np.append(Ksp3, Ksp4, axis=1)

Ksp = np.concatenate((Ksp12, Ksp34), axis=0)


Kps1 = K[2:23, 0:2]

Kps2 = K[24:46, 0:2]

Kps3 = K[2:23, 23]

Kps4 = K[24:46, 23]

Kps3 = Ksp3.reshape((21,1))

Kps4 = Ksp4.reshape((22,1))

Kps12 = np.append(Kps1, Kps2, axis=0)

Kps34 = np.append(Kps3, Kps4, axis=0)

Kps = np.concatenate((Kps12, Kps34), axis=1)


Fp = np.concatenate((F[2:23], F[24:46]))

#calculate displacement and force vectors

iKpp = np.linalg.inv(Kpp)

up = iKpp @ Fp

Fs = Ksp @ up

z1 = np.zeros([2,1])

up = up.reshape((43,1))

U1 = np.append(z1, up[0:21])

z2 = np.zeros([1])

U2 = np.append(U1, z2)

U = np.append(U2, up[21:44])
    


#calculate internal force at element 9

ed = [U[17], U[18], U[19], U[20]]

ed = np.array(ed)

disp = beta[8] @ ed

P = ke_local[8] @ disp

print(P)









    

    