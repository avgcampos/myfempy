# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 21:57:56 2021

PROGRAMA PARA VALIDACAO NUMERICA 
VIGA BI-ENGASTADA COM CARGA DISTRIBUIDA NO COMPRIMENTO

@author: viniv
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

b = 50
h = 100
w = -10
L = 1000
I = (1/12)*b*h**3
E = 200E3


x = np.linspace(0,L,1000)
y = np.linspace(-h/2,h/2,1000)

V = (-w/2)*(L-2*x)

M = (-w/12)*(6*L*x - 6*x**2 - L**2)

v = -((w*x**2)/(24*E*I))*(L-x)**2


plt.close('all')
plt.figure()
plt.subplot(311)
plt.plot(x,V,'r')
plt.grid('on')
plt.subplot(312)
plt.plot(x,M,'b')
plt.grid('on')
plt.subplot(313)
plt.plot(x,v,'k')
plt.grid('on')
plt.show()


M_xy = M*np.ones((len(x),len(y)))
X, Y = np.meshgrid(x, y)
Sig_xx = -(Y/I)*M_xy


plt.figure()
CS = plt.contourf(X, Y, Sig_xx, levels=20, cmap=cm.coolwarm)
plt.contour(CS,colors='k')
plt.xlabel('L')
plt.ylabel('h')
plt.colorbar(CS)







