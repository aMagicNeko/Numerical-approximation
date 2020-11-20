#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:08:55 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

import numpy
import math
from scipy import optimize
from sympy import symbols
from sympy import sympify
from sympy import exp
from sympy import integrate
from sympy import simplify
from sympy import lambdify
from scipy import integrate
from sympy import pi
from sympy import sqrt
import matplotlib.pyplot as plt
def chebyshev(n):
    """
    Create Chebyshev Orthogonal Polynomials

    Parameters
    ----------
    n : int
        the highest order

    Returns
    -------
    P : list
        Chebyshev polynomials

    """
    x=symbols('x',real=True)
    P=[sympify(1),sympify(x)]
    for k in range(2,n+1):
        P.append(2*x*P[k-1]-P[k-2])
    return P

def square_approximation(f,n):
    """
    Create best square approximation on [-1,1] by using Chebyshev polynomials

    Parameters
    ----------
    f : sympy function,its symbol is x
        
    n : int
        the highst order

    Returns
    -------
    s : list
        the target polynomials which order from 0 to n

    """
    x=symbols('x',real=True)
    P=chebyshev(n)
    f_Dot_P=[]
    for k in range(0,n+1):
        f_Dot_P.append(integrate.quad(lambdify(x,f*P[k]/sqrt(1-x**2)),-1,1)[0])
    s=[]
    s.append(sympify(1/pi*f_Dot_P[0]*P[0]))
    for k in range(1,n+1):
        s.append(s[k-1]+2/pi*f_Dot_P[k]*P[k])
    return s

if __name__=="__main__":
    #f=1/(1+25*x**2) on [0,1]
    x=symbols('x',real=True)
    f=1/(1+25*x**2)
    f=f.subs(x,(x+1)/2) #Coordinate transformation
    s=square_approximation(f,11)
    temp=[exp.subs(x,2*x-1) for exp in s]
    s=temp
    xx=numpy.linspace(0,1,1000)
    realy=lambdify(x,f.subs(x,2*x-1))(xx)
    y1=lambdify(x,s[5])(xx)
    plt.figure(1)
    plt.plot(xx,numpy.frompyfunc(abs,1,1)(y1-realy),label="k=5")
    plt.xlabel("x")
    plt.ylabel("$\delta_{k}$")
    plt.title("Chebyshev")
    plt.legend()
    y2=lambdify(x,s[10])(xx)
    plt.figure(2)
    plt.plot(xx,numpy.frompyfunc(abs,1,1)(y2-realy),label="k=10")
    plt.xlabel("x")
    plt.ylabel("$\delta_{k}$")
    plt.title("Chebyshev")
    plt.legend()