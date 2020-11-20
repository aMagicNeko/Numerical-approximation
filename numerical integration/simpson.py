#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 12:58:18 2020

@author: aMagicNeko

email: smz129@outlook.com
"""

import numpy as np
import math
from scipy import integrate
def simpson(f,a,b):
    """
    Compute the itegration of f from a to b

    Parameters
    ----------
    f : function
    
    a : start of the interval
    
    b : end of the interval
    
    n : number of segments
    
    return
    ----------
    value of the integration
    """
    return (b-a)/6*(f(a)+4*f((a+b)/2)+f(b))

def piecewise_simpson(f,a,b,n):
    """
    use piecewise simpson method to compute the itegration of f from a to b
    
    Parameters
    ----------
    f : function
    
    a : start of the interval
    
    b : end of the interval
    
    n : number of segments
    
    return
    ----------
    value of the integration
    """
    x=[]
    h=(b-a)/n
    for k in range(0,n+1):
        x.append(a+h*k)
    res=2*sum([f(y) for y in x[1:n]])
    res+=f(a)+f(b)
    for k in range(1,n+1):
        res+=4*f((x[k-1]+x[k])/2)
    res*=h/6
    return res

def piecewiseTrapezoid(f,a,b,n):
    """
    

    Parameters
    ----------
    f : function
    
    a : start of the interval
    
    b : end of the interval
    
    n : number of segments
    
    return
    ----------
    value of the integration

    """
    h=(b-a)/n
    temp=2*sum([f(a+i*h) for i in range(1,n)])
    res=(temp+f(a)+f(b))*h/2
    return res
if __name__=="__main__":
    def N(x):
        M=200
        miu=1.7
        sigma=0.1
        return M/(sigma*math.sqrt(2*math.pi))*math.exp(((x-miu)**2)/(-2*(sigma**2)))
    a=1.8
    b=1.9
    res1=simpson(N,a,b)
    res2=piecewise_simpson(N,a,b,10)
    res=integrate.quad(N,a,b)
    Er1=abs(res1-res[0])/res[0]
    Er2=abs(res2-res[0])/res[0]
