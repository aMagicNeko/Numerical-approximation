#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:29:30 2020

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
def inner_product(xx,ff,gg):
    L=[ff(y)*gg(y) for y in xx]
    return sum(L)
def inner_product_with_f(xx,fx,ff):
    N=len(xx)
    L=[ff(xx[y])*fx[y] for y in range(0,N)]
    return sum(L)

def squares_least(xx,fx,n):
    """
    Create squares least approximation polynomials

    Parameters
    ----------
    xx : ndarray shape (N)
        x of data points
    fx : ndarray shape (N)
        y of data points 
    n : int
        highest order of target polynomials

    Returns
    -------
    Fit : list
        target polynomials, the index is its order

    """
    N=len(xx)
    #Create orthogonal polynomials
    phi=[1,x-sum(xx)/N]
    Phi=[lambdify(x,y) for y in phi]
    xPhi=[lambdify(x,y*x) for y in phi]
    norm_phi=[]
    norm_phi.append(inner_product(xx,Phi[0],Phi[0]))
    norm_phi.append(inner_product(xx,Phi[1],Phi[1]))
    for k in range(2,n+1):
        alph=inner_product(xx,xPhi[k-1],Phi[k-1])/norm_phi[k-1]
        beta=norm_phi[k-1]/norm_phi[k-2]
        phi.append((x-alph)*phi[k-1]-beta*phi[k-2])
        Phi.append(lambdify(x,phi[k]))
        xPhi.append(lambdify(x,x*phi[k]))
        norm_phi.append(inner_product(xx,Phi[k],Phi[k]))
    #Compute target polynomials
    alph=[inner_product_with_f(xx,fx,Phi[k])/norm_phi[k] for k in range(0,n+1)]
    Fit=[]
    temp=0
    for k in range(0,n+1):
        temp+=alph[k]*phi[k]
        Fit.append(temp)
    return Fit





if __name__=="__main__":
    import matplotlib.pyplot as plt
    x=symbols('x',real=True)
    xx=numpy.arange(5,45,5)
    delta=[3.42,5.96,31.14,41.76,74.54,94.32,133.78,169.16]
    Fit=squares_least(xx,delta,10)
    plt.figure(1)
    xxx=numpy.linspace(0,40,1000)
    for k in range(1,7):
        plt.subplot(3,2,k)
        fit=numpy.frompyfunc(lambdify(x,Fit[k]),1,1)(xxx)
        plt.plot(xxx,fit,label='least squares ')
        plt.plot(xx,delta,"r*",label='real')
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('alpha='+str(k))