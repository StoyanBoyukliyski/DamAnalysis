# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate


#Original points - Keep those to reset
'''
y1 = 0
y2 = 2
y3 = 4 
y4 = 6
y5 = 8 #Keep this constant
x1 = 0
x2 = 2
x3 = 5
x4 = 8.5
x5 = 10 #Keep this constant
'''

#Points that define the curvature of the dam (Feel free to play around)
y1 = 0
y2 = 2
y3 = 4 
y4 = 6
y5 = 8 #Keep this constant
x1 = 0
x2 = 2
x3 = 5
x4 = 8.5
x5 = 10 #Keep this constant

#Water level (Change to see difference)
waterlevel = 3

#Properties of the water
rho = 1000
g = 9.81

#Properties of the Dam
Ec = 25*10**9
vc = 0.22

#Vectorization in calculating the functional coefficients 
Y = np.array([y1,y2,y3,y4,y5])
X = np.array([[1,x1,x1**2,x1**3,x1**4],
              [1,x2,x2**2,x2**3,x2**4],
              [1,x3,x3**2,x3**3,x3**4],
              [1,x4,x4**2,x4**3,x4**4],
              [1,x5,x5**2,x5**3,x5**4]])


#Finding the Coefficients of the curvature equation
A = np.matmul(np.linalg.inv(X),Y)




def InternalStresses(x,z,t=0.1):
    #Calculate internal stresses at each cross-section
    y = (A[0]+A[1]*x+A[2]*x**2+A[3]*x**3+A[4]*x**4)
    ynew = 8 - y
    I = t*(ynew)**3/12
    w1 = rho*g*waterlevel*t
    w2 = 0
    x1c = 0
    x2c = waterlevel
    
    #Vectorize solution
    W = np.array([w1,w2])
    X = np.array([[1,x1c],
                  [1,x2c]])
    C = np.matmul(np.linalg.inv(X),W)
    
    #Calculate based on Bernoulli equation
    if x > 0 and x <waterlevel:
        M = C[0]*x**2/2 + C[1]*x**3/6 - (C[0]*waterlevel +C[1]*waterlevel**2/2)*x - (C[0]*waterlevel**2/2+C[1]*waterlevel**3/6-(C[0]*waterlevel+C[1]*waterlevel**2/2)*waterlevel)
        V = C[0]*x + C[1]*x**2/2 - (C[0]*waterlevel+ C[1]*waterlevel**2/2)
    else:
        V = 0
        M = 0
        
    #Estimate the centroid for sanity check
    centroid = 4 +y/2
    
    #Calculate the first moment of area
    vec1 = np.array([0, ynew, 0])
    mat1 = np.array([[1, 0, 0],
            [1, ynew/2, ynew**2/8],
            [1, ynew, ynew**2]])
    coef = np.matmul(np.linalg.inv(mat1),vec1)
    Q = t*(coef[0]+coef[1]*(z-y)+coef[2]*(z-y)**2)/2
    
    #Calculate internal shear and normal stresses
    if z > y and z < 8:
        sigma = M* (z-centroid)/ I
        tau = (V*Q)/(I*t)
    else:
        sigma = None
        tau = None
    
    return [sigma,M,tau,V]

#Set the numerical accuracy of the plot
accuracy = 300
u = np.linspace(0,8,accuracy)
v = np.linspace(0,10,accuracy)

#Initiate parameters and iterate stresses in form for plotting
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
        Sigma = InternalStresses(i,j)
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
    
#Plot figures for internal normal stresses perpendicular to dam's cross-section
figure = plt.figure(1)
ax2 = figure.add_subplot(111)
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
ax2.plot([8,10], [waterlevel,waterlevel], "b-")
CS = ax2.contourf(U,V, matrix, levels = contours, extend = "both", cmap = "inferno")
figure.colorbar(CS)
ax2.axis("equal")

#Plot figures for internal shear stresses parallel and perpendicular to dam's cross-section
figure2 = plt.figure(2)
ax2 = figure2.add_subplot(111)

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
ax2.plot([8,10], [waterlevel,waterlevel], "b-")
CS1 = ax2.contourf(U,V, matrixtau, 
                   levels = contours2, extend = "both", cmap = "inferno")
figure2.colorbar(CS1)
ax2.axis("equal")

