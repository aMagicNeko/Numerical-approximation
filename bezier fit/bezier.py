#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 13:22:31 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d
from scipy.special import comb, perm
def bezier(x,t):
    """
    Use de Casteljau algorithm

    Parameters
    ----------
    x : numpy.array
        data points
    t : float
        point to get

    Returns
    -------
    xx : points represented by t

    """
    P=x.copy()
    k=len(x)-1
    while k >= 0:
        k-=1
        for j in range(0,k+1):
            P[j]=(1-t)*P[j]+t*P[j+1]
    return P[0]
def div(x,t):
    """
    Divide the Bezier curve into two pieces

    Parameters
    ----------
    x : ndarray
        data points
    t : float

    Returns
    -------
    ndarray,ndarray
    the control points
    """
    P=x.copy()
    k=len(x)-1
    res1=[x[0]]
    res2=[x[-1]]
    while k >= 0:
        k-=1
        for j in range(0,k+1):
            P[j]=(1-t)*P[j]+t*P[j+1]
        res1.append(P[0].copy())
    return np.array(res1),P
def ascend_order(x,m):
    """
    

    Parameters
    ----------
    xy : ndarray
        data points
    m : int
        increase order

    Returns
    -------
    nx : new points 
        ndarray

    """
    n=len(x)-1 #the present order
    nx=np.zeros((n+m+1,2))
    for i in range(0,n+m+1):
        for j in range(max(i-m,0),min(i+1,n+1)):
            nx[i]+=comb(n,j)*comb(m,i-j)/comb(n+m,i)*x[j]
    return nx






if __name__=="__main__":
    x=np.array([1,4,5,6])
    y=[2,4,8,3]
    plt.plot(x,y,'b')
    x1=x+1
    y1=[13,14,16,14]
    plt.plot(x1,y1,'b')
    plt.plot([1,2],[2,13],'b')
    plt.plot([6,7],[3,14],'b')
    x=np.array(list(x)+list(x1))
    y=np.array(y+y1)
    xy=np.zeros((8,2))
    xy[:,0]=x
    xy[:,1]=y
    T=np.linspace(0,1,10000)
    res=[]
    for t in T:
        res.append(bezier(xy,t))
    res=np.array(res)
    plt.plot(res[:,0],res[:,1])
    t0=0.14
    p0=bezier(xy,t0)
    plt.plot(p0[0],p0[1],'b*')
    dp1,dp2=div(xy,t0)
    plt.plot(dp1[:-2,0],dp1[:-2,1],"r+")
    plt.plot(dp2[1:,0],dp2[1:,1],"o")
    plt.plot(dp1[:,0],dp1[:,1],"--",color="deeppink")
    plt.plot(dp2[:,0],dp2[:,1],"--",color="deeppink")
    nx=ascend_order(xy,2)
    plt.plot(nx[:,0],nx[:,1],'v',color="darkgray")