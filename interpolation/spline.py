#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 23:08:44 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

import numpy
import matplotlib.pyplot as plt
def spline(x,y,y1,y2,xx):
    """
    spline interpolation

    Parameters
    ----------
    x : ndarray shape: (n)
        x of data points
    y : ndarray shape: (n)
        y of data points
    y1 : float
        first derivation of f at x[0]
    y2 : float
        first derivation of f at x[-1]
    xx : ndarray shape: (n)
        x of target points

    Returns
    -------
    yy : ndarray shape:(n)
        y of target points

    """
    yy=numpy.zeros(len(xx))
    n=len(x)
    h=numpy.zeros(n-1)
    for k in range(0,n-1):
        h[k]=x[k+1]-x[k]
    Df=numpy.zeros((n,3))
    Df[:,0]=y
    for k in range(1,3):
        for j in range(k,3):
            Df[j,k]=(Df[j-1,k-1]-Df[j,k-1])/(x[j-k]-x[j])
    d=numpy.zeros(n)
    d[0]=6/h[0]*(Df[1,1]-y1)
    d[n-1]=6/h[n-2]*(y2-Df[n-1,1])
    for k in range(1,n-1):
        d[k]=6*Df[k+1,2]
    A=numpy.zeros((n,n))
    A[0,0]=2
    A[0,1]=1
    A[n-1,n-1]=2
    A[n-1,n-2]=1
    for k in range(1,n-1):
        A[k,k]=2
        A[k,k-1]=h[k-1]/(h[k]+h[k-1])
        A[k,k+1]=h[k]/(h[k]+h[k-1])
    M=numpy.linalg.solve(A,d)
    for k in range(0,n-1):
        flag=((xx>=x[k]) & (xx<=x[k+1]))
        yy[flag]=((x[k+1]-xx[flag])**3)/(6*h[k])*M[k]+((xx[flag]-x[k])**3)/(6*h[k])*M[k+1]\
            +(y[k]-M[k]/6*(h[k]**2))*(x[k+1]-xx[flag])/h[k]\
                +(y[k+1]-M[k+1]/6*(h[k]**2))*(xx[flag]-x[k])/h[k]
    return yy



if __name__=="__main__":
    N=numpy.array([4,10,14,20])
    xx=numpy.linspace(-5,5,1000)
    realy=1/(1+xx**2)
    k=1
    for n in N:
        step=10/n
        x=numpy.arange(-5,5+step,step)
        y=1/(1+x**2)
        y1=((-2)*x[0])/((1+x[0]**2)**2)
        y2=((-2)*x[-1])/((1+x[-1]**2)**2)
        yy=spline(x,y,y1,y2,xx)
        plt.figure(1)
        plt.subplot(2,2,k)
        k=k+1
        plt.plot(xx,yy,"r",label="Spline")
        plt.plot(xx,realy,"b--",label="real")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.title("n="+str(n))