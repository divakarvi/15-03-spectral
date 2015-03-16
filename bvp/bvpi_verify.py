#! /usr/bin/env python
import numpy as np
from scipy import linalg
from StringIO import StringIO
from pprint import pprint
import sys
import os

def read(fname):
    return np.genfromtxt(fname, dtype='float')

def intcheb(f):
    n = len(f)
    ff = np.zeros(n)
    for i in range(1,n-1):
        ff[i] = (f[i-1]-f[i+1])/(2*i)
    return ff
    
def error(a, u, du, f):
    e = du - a*a*intcheb(u) - intcheb(f)
    e[0] = 0;
    #print 'du = ', du
    #print 'u = ', u
    #print 'int of u = ', intcheb(u)
    #print 'rem = ', du -a*a*intcheb(u)
    if sum(abs(f))> 1e-10:
        return linalg.norm(e)/linalg.norm(f)
    else:
        return linalg.norm(e)


if(__name__=="__main__"):
    fp = open("DBG/a.txt", 'r');
    a = float(fp.read())
    u = read("DBG/u.txt")
    du = read("DBG/du.txt")
    f = read("DBG/f.txt")
    e = error(a, u, du, f)
    print '\t\t error = ', e
