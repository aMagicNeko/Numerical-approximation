#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:04:10 2020

@author: aMagicNeko

email: smz129@outlook.com
"""
#use x^i as basis function to do best square approximation
import numpy
import math
from sympy import *
import scipy.linalg
import matplotlib.pyplot as plt
x=symbols('x',real=True)
f=1/(1+25*x**2)
P=[1,x]
for k in range(2,11):
    P.append(x*P[k-1])
H1=scipy.linalg.hilbert(6)
H2=scipy.linalg.hilbert(11)
d1=numpy.zeros(6)
d2=numpy.zeros(11)
for k in range(0,11):
    d2[k]=integrate(f*P[k],\
                    (x,0,1)).evalf()
d1=d2[:6]
A1=scipy.linalg.solve(H1,d1)
A2=scipy.linalg.solve(H2,d2)
Fit1=0
Fit2=0
xx=numpy.linspace(0,1,1000)
for k in range(0,11):
    Fit2+=A2[k]*P[k]
for k in range(0,6):
    Fit1+=A1[k]*P[k]
fit2=lambdify(x,Fit2)
fit1=lambdify(x,Fit1)
realy=lambdify(x,f)(xx)
y1=fit1(xx)
y2=fit2(xx)
plt.figure(1)
plt.plot(xx,numpy.frompyfunc(abs,1,1)\
         (y1-realy),label="k=5")
plt.xlabel("x")
plt.ylabel("$\delta_{k}$")
plt.legend()
plt.figure(2)
plt.plot(xx,numpy.frompyfunc(abs,1,1)\
         (y2-realy),label="k=10")
plt.legend()
plt.xlabel("x")
plt.ylabel("$\delta_{k}$")
print(Fit1)
print(Fit2)