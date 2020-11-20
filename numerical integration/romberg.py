#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 13:04:11 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

import numpy as np
import math
from scipy import integrate
def romberg(f,a,b,delta):
    """
    Use romberg method to compute the integration of f from a to b

    Parameters
    ----------
    f : funcion
    
    a : float
        start of the interval
    b : float
        end of the interval
    delta :float
        bound of the error

    Returns
    -------
    float
        value of integration.

    """
    h=b-a
    k=0
    T=np.zeros((2,10000000))
    T[0,0]=h/2*(f(a)+f(b))
    while True:
        k+=1
        h/=2
        index=k%2
        temp=sum([f(a+i*h) for i in range(1,2**k)])
        T[index,0]=h/2*(f(a)+f(b)+2*temp)
        for j in range(1,k+1):
            T[index,j]=1/(4**j-1)*(4**j*T[index,j-1]-T[(k+1)%2,j-1])
        err=abs(T[index,k]-T[k%2,k-1])
        if err<delta:
            return T[index,k]


if __name__=="__main__":
    def f(x):
        return math.log(x)
    res1=romberg(f,1,2,1e-6)
    res=integrate.quad(f,1,2)
    err=abs((res[0]-res1)/res[0])