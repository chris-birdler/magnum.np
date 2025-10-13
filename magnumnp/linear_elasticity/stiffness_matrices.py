#
# This file is part of the magnum.np distribution
# (https://gitlab.com/magnum.np/magnum.np).
# Copyright (c) 2023 magnum.np team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np

__all__ = ["C_cubic", "C_isotropic", "C_trigonal", "RotateC"]

def C_cubic(C11, C12, C44):
    C = [[C11, C12, C12,   0,   0,   0],
         [C12, C11, C12,   0,   0,   0],
         [C12, C12, C11,   0,   0,   0],
         [  0,   0,   0, C44,   0,   0],
         [  0,   0,   0,   0, C44,   0],
         [  0,   0,   0,   0,   0, C44]]

    return np.array(C)
    
def C_isotropic(E, nu):
    C11 = E*(1.-nu)/((1.+nu)*(1.-2.*nu))
    C12 = E*nu/((1.+nu)*(1.-2.*nu))
    C44 = 0.5*(C11 - C12)

    return C_cubic(C11, C12, C44)

def C_trigonal(C11, C12, C13, C14, C33, C44):
    C = [[C11,  C12,  C13,  C14,   0,   0],
         [C12,  C11,  C13, -C14,   0,   0],
         [C13,  C13,  C33,    0,   0,   0],
         [C14, -C14,    0,  C44,   0,   0],
         [  0,    0,    0,    0, C44, C14],
         [  0,    0,    0,    0, C14, 0.5*(C11-C12)]]

    return np.array(C)

def RotateC(C, axis, angle):
    R = RMat(axis, angle)

    map = np.array([[0,5,4],[5,1,3],[4,3,2]])

    Ct = np.zeros((3,3,3,3))
    Cr = np.zeros((3,3,3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Ct[i,j,k,l] = C[map[i,j],map[k,l]]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Cr[i,j,k,l] += np.dot(R[l], Ct[i,j,k,:])
    Ct = Cr
    Cr = np.zeros((3,3,3,3))

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Cr[i,j,k,l] += np.dot(R[k], Ct[i,j,:,l])
    Ct = Cr
    Cr = np.zeros((3,3,3,3))

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Cr[i,j,k,l] += np.dot(R[j], Ct[i,:,k,l])
    Ct = Cr 
    Cr = np.zeros((3,3,3,3))

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Cr[i,j,k,l] += np.dot(R[i], Ct[:,j,k,l])
    
    res = np.zeros((6,6))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    res[map[i,j],map[k,l]] = Cr[i,j,k,l]

    return res

def cross_product_matrix(u):
    ucm = np.zeros((3,3))
    ucm[0,1] = -u[2]
    ucm[0,2] = u[1]
    
    ucm[1,0] = u[2]
    ucm[1,2] = -u[0]
    
    ucm[2,0] = -u[1]
    ucm[2,1] = u[0]

    return ucm

def RMat(axis, angle):
    if axis == "x":
        u = np.array([1,0,0])
    elif axis == "y":
        u = np.array([0,1,0])
    elif axis == "z":
        u = np.array([0,0,1])
    else:
        u = axis

    ucm = cross_product_matrix(u)
    uom = np.outer(u,u)

    R = np.cos(angle)*np.eye(3) 
    R += np.sin(angle)*ucm 
    R += (1.-np.cos(angle))*uom

    return R