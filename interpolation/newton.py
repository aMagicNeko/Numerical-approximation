#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 22:09:42 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

import numpy
import math
import matplotlib.pyplot as plt


def newton(x,y,xx):
    """
    Newton interpolation

    Parameters
    ----------
    x : ndarray shape:(n)
        x of data points
    y : ndarray shape:(n)
        y of data points
    xx : ndarray shape:(k)
        x of target points

    Returns
    -------
    yy : ndarray shape:(k)
        y of target points

    """
    n=len(x)
    Df=numpy.zeros((n,n))
    Df[:,0]=y
    for k in range(1,n):
        for j in range(k,n):
            Df[(j,k)]=(Df[(j-1,k-1)]-Df[(j,k-1)])/(x[j-k]-x[j])
    yy=numpy.zeros(len(xx))
    products=[]
    products.append(numpy.zeros(len(xx))+1)
    for k in range(1,n):
        products.append(products[k-1]*(xx-x[k-1]))
    products=numpy.array(products)
    for k in range(0,n):
        yy+=products[k,:]*Df[k,k]
    return yy


if __name__=="__main__":
    N=numpy.array([4,10,14,20])
    xx=numpy.linspace(-5,5,1000)
    realy=1/(1+xx**2)
    #use equidistant nodes to interpolate
    k=1
    for n in N:
        step=10/n
        x=numpy.arange(-5,5+step,step)
        y=1/(1+x**2)
        yy=newton(x,y,xx)
        plt.figure(1)
        plt.subplot(2,2,k)
        k=k+1
        plt.plot(xx,yy,"r",label="newton")
        plt.plot(xx,realy,"b--",label="real")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.title("n="+str(n))
    #use Chebychev points to interpolate
    k=1
    for n in N:
        temp=numpy.array([math.pi/2*(2*i+1)/(n+1) for i in range(0,n+1)])
        z=numpy.cos(temp)
        x=10/2*z
        y=1/(1+x**2)
        yy=newton(x,y,xx)
        plt.figure(2)
        plt.subplot(2,2,k)
        k=k+1
        plt.plot(xx,yy,"r",label="newton")
        plt.plot(xx,realy,"b--",label="real")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.title("n="+str(n))
        