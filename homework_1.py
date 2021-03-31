#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 01:56:44 2021

@author: Tomy
"""

from numpy import array , pi , zeros , ix_ , sin , cos
from beam_element import beam_element
import matplotlib.pylab as plt

densidad_hormigon = 2500 #kg/m3

g = 9.8



xy = array ([[0,0],
             [0,3],
             [0,5],
             [0,6],
             [6,6.5],
             [6,5],
             [6,3],
             [6,0]
             ])

conec = array ([[0,1],
                [1,2],
                [2,3],
                [7,6],
                [6,5],
                [5,4],
                [1,6],
                [2,5],
                [3,4]
                ],
                dtype = int
                )


#COLUMNS. Elements: 0,1,2,3,4,5
b0 = 30e-2
h0 = 30e-2
properties_0 = {}
properties_0["E"] = 25.6e9
properties_0["A"] = b0*h0
properties_0["I"] = (b0 * h0**3)/12
properties_0["qx"] = densidad_hormigon * b0*h0 * g
properties_0["qy"] = densidad_hormigon * b0*h0 * g

#BEAMS. Elements: 6,7
b1 = 20e-2
h1 = 40e-2
properties_1 = {}
properties_1["E"] = 25.6e9
properties_1["A"] = b1*h1
properties_1["I"] = (b1 * h1**3)/12
properties_1["qx"] = densidad_hormigon * b1*h1 * g
properties_1["qy"] = densidad_hormigon * b1*h1 * g

#ROOF. Element: 8
b2 = 20e-2
h2 = 20e-2
properties_2 = {}
properties_2["E"] = 25.6e9
properties_2["A"] = b2*h2
properties_2["I"] = (b2 * h2**3)/12
properties_2["qx"] = densidad_hormigon * b2*h2 * g
properties_2["qy"] = densidad_hormigon * b2*h2 * g





properties = [ properties_0 , properties_0 , properties_0 , properties_0 , properties_0 , properties_0 , #COLUMNS
               properties_1 , properties_1 ,                                                             #BEAMS
               properties_2                                                                              #ROOF
               ]

Nnodes = xy.shape[0]
Nelems = conec.shape[0]

Ndofs_per_node = 3
Ndofs = Ndofs_per_node * Nnodes
K = zeros((Ndofs,Ndofs))
f = zeros((Ndofs,1))

for e in range(Nelems):
    ni = conec[e,0]
    nj = conec[e,1]
    
    xy_e = xy[[ ni , nj ],:] #reescribiendo xy
    
    ke , fe = beam_element( xy_e , properties[e])
    
    d = [ 3*ni , 3*ni+1 , 3*ni+2 , 3*nj , 3*nj+1 ,3*nj+2 ]
    
    #DIRECT STIFFNESS METHOD
    for i in range(2*Ndofs_per_node):
        p = d[i]
        for j in range(2*Ndofs_per_node):
            q = d[j]
            K[p,q] += ke[i,j]
        f[p] += fe[i]
        
#HAND CALCULATION OF f

print (f"K={K}")
print (f"f={f}")

plt.figure()
plt.matshow(K)
plt.colorbar()
plt.show()

# SYSTEM PARTITIONING AND SOLUTION

free_DOFs = [3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20]
c_DOFs = [0,1,2,21,22,23]

Kff = K[ix_(free_DOFs, free_DOFs)]
Kfc = K[ix_(free_DOFs, c_DOFs)]
Kcf = K[ix_(c_DOFs, free_DOFs)]
Kcc = K[ix_(c_DOFs, c_DOFs)]

ff = f[free_DOFs]
fc = f[c_DOFs]

#SOLVE
from scipy.linalg import solve
u = zeros((Ndofs,1))

u[free_DOFs] = solve( Kff , ff )

#GET REACTION FORCES
R = Kcf @ u[free_DOFs] + Kcc @ u[c_DOFs] - fc 

print (f"u={u}")
print (f"R={R}")
    
    
    