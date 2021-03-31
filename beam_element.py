#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 00:34:07 2021

@author: Tomy
"""

from sympy import pretty
from numpy import array, zeros, ix_
from scipy.linalg import norm

def beam_element(xy, properties):
    E = properties["E"]
    I = properties["I"]
    A = properties["A"]
    
    qx = properties["qx"]
    qy = properties["qy"]
    
    xi = xy[0,:]
    xj = xy[1,:]
    
    print (f"xi: {xi}")
    print (f"xj: {xj}")
    
    L = norm(xj - xi)
    cos = (xj[0] - xi[0])/L
    sen = (xj[1] - xi[1])/L
    
    ke = zeros((6,6))
    fe = zeros((6,1))
    ke_tilde = zeros((6,6))
    
    ke_tilde[0,0] = A*E / L
    ke_tilde[3,3] = A*E / L
    ke_tilde[0,3] = -A*E / L
    ke_tilde[3,0] = -A*E / L
    
    bending_dofs = ix_([1,2,4,5],[1,2,4,5])
    
    ke_tilde[ bending_dofs ] = E*I*array([[12/L**3, 6/L**2, -12/L**3, 6/L**2],
                                          [6/L**2, 4/L, -6/L**2, 2/L],
                                          [-12/L**3, -6/L**2, 12/L**3, -6/L**2], 
                                          [6/L**2, 2/L, -6/L**2, 4/L]])
    
    T = zeros((6,6))
    T[0:2,0:2] = array ([[cos , -sen ], [sen , cos]])
    T[3:5,3:5] = array ([[cos , -sen ], [sen , cos]])
    T[2,2] = 1.0
    T[5,5] = 1.0
    
    #COMPUTE fe FROM qx AND qy...
    
    if cos == 0 and sen == 1 :
        fe[0] = -qx*L/2
        fe[1] = 0
        fe[2] = 0
        fe[3] = -qx*L/2
        fe[4] = 0
        fe[5] = 0
    else:
        fe[0] = -qy*L/2 * -sen
        fe[1] = -qy*L/2 *  cos
        fe[2] = -qy*L**2/12 
        fe[3] = -qy*L/2 * -sen
        fe[4] = -qy*L/2 *  cos
        fe[5] = +qy*L**2/12 
    
    
    
    
    ke = T @ ke_tilde @ T.T
    
    print (pretty(fe))
    
    return ke, fe
"""
xy = array([[0,0] , 
            [0,1]])

densidad_hormigon = 2500 #kg/m3

g = 9.8

b0 = 30e-2
h0 = 30e-2

properties = {}
properties["E"] = 25.6e9
properties["I"] = (b0 * h0**3)/12
properties["A"] = b0*h0
properties["qx"] = densidad_hormigon * g * b0*h0
properties["qy"] = densidad_hormigon * g * b0*h0


beam_element( xy , properties )
"""

