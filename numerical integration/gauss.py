#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 13:17:02 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

"""

I have not typed the codes of piecewise method of gauss_legendre integration.

"""
from scipy import integrate
import math
import numpy as np
import matplotlib.pyplot as plt
def gauss_legendre(f,n):
    """
    Compute the integration of f on [-1,1]

    Parameters
    ----------
    f : function
        
    n : int
        order of the numerical integration

    Returns
    -------
    res : float
        the result of the integration
    err : float
        relative error

    """
    x,weight=np.polynomial.legendre.leggauss(n)
    res=sum([weight[k]*g(x[k]) for k in range(n)])
    real=integrate.quad(f,-1,1)[0]
    err=abs((res-real)/real)
    return res,err
if __name__=="__main__":
    def f(x):
        return math.sqrt(1+math.exp(x))
    res=res=integrate.quad(f,0,4)
    def g(x):
        return 2*math.sqrt(1+math.exp(2*x+2))
    err3=gauss_legendre(g,3)[1]
    err5=gauss_legendre(g,5)[1]
    err7=gauss_legendre(g,7)[1]
    plt.plot([3,5,7],[err3,err5,err7])
