#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 22:59:43 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

import numpy
import matplotlib.pyplot as plt
import math

def _hermite(x,y,y1,xx):
    res=(1-2*(xx-x[0])/(x[0]-x[1]))*(((xx-x[1])/(x[0]-x[1]))**2)*y[0]+\
    (xx-x[0])*(((xx-x[1])/(x[0]-x[1]))**2)*y1[0]+\
    (1-2*((xx-x[1])/(x[1]-x[0])))*(((xx-x[0])/(x[1]-x[0]))**2)*y[1]+\
    (xx-x[1])*(((xx-x[0])/(x[1]-x[0]))**2)*y1[1]
    return res
    
def hermite(x,y,y1,xx):
    """
    Hermite interpolation.

    Parameters
    ----------
    x : ndarray shape: (n)
        x of data points
    y : ndarray shape: (n)
        y of data points
    y1 : ndarray shape: (n)
        first derivation of f(x)
    xx : ndarray shape: (k)
        x of target points
    Returns
    -------
    yy : ndarray shape: (k)
        y of target points 

    """
    yy=numpy.zeros(len(xx))
    n=len(x)
    for k in range(0,n-1):
        flag=((xx>=x[k]) & (xx<=x[k+1]))
        yy[flag]=_hermite(x[k:k+2],y[k:k+2],y1[k:k+2],xx[flag])
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
        y1=((-2)*x)/((1+x**2)**2)
        yy=hermite(x,y,y1,xx)
        plt.figure(1)
        plt.subplot(2,2,k)
        k=k+1
        plt.plot(xx,yy,"r",label="Hermite")
        plt.plot(xx,realy,"b--",label="real")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.title("n="+str(n))