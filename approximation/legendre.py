#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:09:53 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

"""

Use Legendre Orthogonal Polynomials to do Best square approximation

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
import matplotlib.pyplot as plt
def legendre(n):
    """
    Create Legendre Orthogonal Polynomials

    Parameters
    ----------
    n : int
        the highest order

    Returns
    -------
    P : list
        Legendre polynomials

    """
    x=symbols('x',real=True)
    P=[sympify(1),x]
    for k in range(2,n+1):
        P.append(((2*k-1)*x*P[k-1]-(k-1)*P[k-2])/k)
    for i in range(0,n+1):
        P[i]=simplify(P[i])
    return P
def square_approximation(f,n):
    """
    Create best square approximation on [-1,1] by using Legendre polynomials

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
    P=legendre(n)
    f_Dot_P=[]
    for k in range(0,n+1):
        f_Dot_P.append(integrate(f*P[k],(x,-1,1)))
    s=[]
    s.append(sympify(1/2*f_Dot_P[0]*P[0]))
    for k in range(1,n+1):
        s.append(s[k-1]+(2*k+1)/2*f_Dot_P[k]*P[k])
    return s
if __name__=="__main__":
    #approximate e^x at order 1,3,5,10
    from matplotlib import pyplot as plt
    x=symbols("x",real=True)
    f=exp(x)
    s=square_approximation(f,10)
    N=[1,3,5,10]
    xx=numpy.linspace(-1,1,1000)
    i=1
    plt.figure(1)
    f_l=lambdify(x,f)
    for k in N:
        g=lambdify(x,s[k])
        yy=numpy.array([abs(g(y)-f_l(y)) for y in xx ])
        plt.subplot(2,2,i)
        i+=1
        plt.plot(xx,yy)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(k)
    #Plot maximum of the error with its order 
    from scipy import optimize
    res=[]
    K=numpy.array(list(range(0,11)))
    for k in K:
        def F(xx):
            g=lambdify(x,s[k])
            if(xx>=-1 and xx<=1):
                return -abs(g(xx)-f_l(xx))
            else:
                return 0
        res.append(optimize.minimize(F,0,options={'disp':False}))
    delta=[]
    for R in res:
        delta.append(-R['fun'])
        #print(-R['fun']) 
    plt.figure(2)
    plt.plot(K,delta)
    plt.xlabel('k')
    plt.ylabel('$\Delta_{k}$')
    
