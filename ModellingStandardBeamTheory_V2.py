# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate

y1 = 0
y2 = 2
y3 = 4 
y4 = 6
y5 = 8
x1 = 0
x2 = 2
x3 = 5
x4 = 8.5
x5 = 10
roh = 1000
g = 9.81
Ec = 25*10**9
vc = 0.22

Y = np.array([y1,y2,y3,y4,y5])
X = np.array([[1,x1,x1**2,x1**3,x1**4],
              [1,x2,x2**2,x2**3,x2**4],
              [1,x3,x3**2,x3**3,x3**4],
              [1,x4,x4**2,x4**3,x4**4],
              [1,x5,x5**2,x5**3,x5**4]])

A = np.matmul(np.linalg.inv(X),Y)

def CreatePlot(x,z,t=0.1):
    y = (A[0]+A[1]*x+A[2]*x**2+A[3]*x**3+A[4]*x**4)
    ynew = 8 - y
    I = t*(ynew)**3/12
    w1 = 78480*t
    w2 = 0
    x1c = 0
    x2c = 8
    
    W = np.array([w1,w2])
    X = np.array([[1,x1c],
                  [1,x2c]])
    C = np.matmul(np.linalg.inv(X),W)
    if x > 0 and x <8:
        M = C[0]*x**2/2 + C[1]*x**3/6 - (C[0]*8 +C[1]*8**2/2)*x -(C[0]*8**2/2+C[1]*8**3/6-(C[0]*8+C[1]*8**2/2)*8)
        V = C[0]*x + C[1]*x**2/2 - (C[0]*8+ C[1]*8**2/2)
    else:
        V = 0
        M = 0
        
    centroid = 4 +y/2
    vec1 = np.array([0, ynew, 0])
    mat1 = np.array([[1, 0, 0],
            [1, ynew/2, ynew**2/8],
            [1, ynew, ynew**2]])
    coef = np.matmul(np.linalg.inv(mat1),vec1)
    Q = t*(coef[0]+coef[1]*(z-y)+coef[2]*(z-y)**2)/2
    if z > y and z < 8:
        sigma = M* (z-centroid)/ I
        tau = (V*Q)/(I*t)
    else:
        sigma = None
        tau = None
    
    return [sigma,M,tau,V]


figure, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize =(15,10))
ax1.plot([A[0] + A[1]*x + A[2]*x**2 +A[3]*x**3 + A[4]*x**4 
          for x in np.arange(0.0,10.0,0.05)],
            [x for x in np.arange(0.0,10.0,0.05)], "r-")
ax1.plot([8,8], [10,0], "r-")
ax1.plot([0,8], [0,0], "r-")
ax1.plot([-2,10], [0,0], "g-")
ax1.plot([8,10], [8,8], "b-")
ax1.axis("equal")

accuracy = 300
u = np.linspace(0,8,accuracy)
v = np.linspace(0,10,accuracy)
vector = []
matrix = []
momentvec= []
momentmatrix = []
vectortau = []
matrixtau = []
sigmamax = 0
taumax = 0
sigmamin = 0
taumin = 0
for i in v:
    for j in u:
        Sigma = CreatePlot(i,j)
        sigma = Sigma[0]
        Moment = Sigma[1]
        newtau = Sigma[2]
        vector.append(sigma)
        momentvec.append(Moment)
        vectortau.append(newtau)
        if sigma != None:
            if sigmamax < sigma:
                sigmamax = sigma
            else:
                pass
            if sigmamin > sigma:
                sigmamin = sigma
            else:
                pass
        if newtau != None:
            if taumax < newtau:
                taumax = newtau
            else:
                pass
        
            if taumin > newtau:
                taumin = newtau
            else:
                pass
        
    momentmatrix.append(momentvec)
    matrix.append(vector)
    matrixtau.append(vectortau)
    vector = []
    vectortau = []
    momentvec = []
    
U,V = np.meshgrid(u,v)
contours = np.linspace(sigmamin, sigmamax ,30)
ax2.plot([A[0] + A[1]*x + A[2]*x**2 +A[3]*x**3 + A[4]*x**4 
          for x in np.arange(0.0,10.0,0.05)],[x
                            for x in np.arange(0.0,10.0,0.05)], "r-")
ax2.plot([(8-A[0] + A[1]*x + A[2]*x**2 +A[3]*x**3 + A[4]*x**4)/2 
          for x in np.arange(0.0,10.0,0.05)],[x 
                            for x in np.arange(0.0,10.0,0.05)], "b--")
ax2.plot([8,8], [10,0], "r-")
ax2.plot([0,8], [0,0], "r-")
ax2.plot([-2,10], [0,0], "g-")
CS = ax2.contourf(U,V, matrix, 
                  levels = contours, extend = "both", cmap = "RdGy")
figure.colorbar(CS)
ax2.axis("equal")


figure2, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize =(20,10))

ax1.plot([A[0] + A[1]*x + A[2]*x**2 +A[3]*x**3 + A[4]*x**4 
          for x in np.arange(0.0,10.0,0.05)],[x 
                            for x in np.arange(0.0,10.0,0.05)], "r-")
ax1.plot([8,8], [10,0], "r-")
ax1.plot([0,8], [0,0], "r-")
ax1.plot([-2,10], [0,0], "g-")
ax1.plot([8,10], [8,8], "b-")
ax1.axis("equal")


contours2 = np.linspace(taumin, taumax,30)
ax2.plot([A[0] + A[1]*x + A[2]*x**2 +A[3]*x**3 + A[4]*x**4 
          for x in np.arange(0.0,10.0,0.05)],[x 
                            for x in np.arange(0.0,10.0,0.05)], "r-")
ax2.plot([(8-A[0] + A[1]*x + A[2]*x**2 +A[3]*x**3 + A[4]*x**4)/2 
          for x in np.arange(0.0,10.0,0.05)],[x 
                            for x in np.arange(0.0,10.0,0.05)], "b--")
ax2.plot([8,8], [10,0], "r-")
ax2.plot([0,8], [0,0], "r-")
ax2.plot([-2,10], [0,0], "g-")
CS1 = ax2.contourf(U,V, matrixtau, 
                   levels = contours2, extend = "both", cmap = "RdGy")
figure2.colorbar(CS1)
ax2.axis("equal")

