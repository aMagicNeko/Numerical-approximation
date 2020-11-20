#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 22:09:42 2020

@author: aMagicNeko

email: smz129@outlook.com
"""
import numpy
import matplotlib.pyplot as plt
import math

def lagrange(x,y,xx):
    """
    Lagrange Interpolation

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
    yy=numpy.zeros(len(xx))
    i=0
    while i<n:
        temp=numpy.zeros(len(xx))+1
        for k in list(range(0,n)):
            if k==i:
                continue
            temp=temp*(xx-x[k])/(x[i]-x[k])
        yy=yy+temp*y[i]
        i=i+1
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
        yy=lagrange(x,y,xx)
        plt.figure(1)
        plt.subplot(2,2,k)
        k=k+1
        plt.plot(xx,yy,"r",label="lagrange")
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
        yy=lagrange(x,y,xx)
        plt.figure(2)
        plt.subplot(2,2,k)
        k+=1
        plt.plot(xx,yy,"r",label="lagrange&chebychev")
        plt.plot(xx,realy,"b--",label="real")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.title("n="+str(n))