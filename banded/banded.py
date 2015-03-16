#! /usr/bin/env python
import numpy as np
from scipy import linalg
from StringIO import StringIO
from pprint import pprint
import sys

def read(fname):
    return np.genfromtxt(fname, dtype='float')

def verify(af, bf, xf, l=1, u=1):
    A = read(af)
    b = read(bf)
    x = read(xf)
    xx = linalg.solve_banded((l,u), A, b)
    #print 'A = \n', A
    #print 'b = \n', b
    #print 'x = \n', x
    #print 'xx = \n', xx
    #print 'shape of x, b: ', x.shape()
    error = linalg.norm(xx-x)/linalg.norm(x)
    print 'rel error = ', error


def test1():
    af = StringIO('0 1 2 \n 3 4 5 \n 6 7 0')
    bf = StringIO('1 1 \n 2 2 \n 3 3');
    xf = StringIO('-0.16666666666 -0.16666666666\n 1.5  1.5\n -1.5 -1.5')
    verify(af, bf, xf, l=1, u=1)

def test2():
    af = StringIO('0 0 1 2 3 \n 0 1 2 3 4 \n 1 2 3 4 5 \n 1 2 3 4 0 \n 1 2 3 0 0');
    bf = StringIO('1 \n 2 \n 3 \n 4 \n 5')
    xf = StringIO('-4 \n -1 \n 6 \n -2 \n -1')
    verify(af, bf, xf, l=2, u=2)

def test3():
    a = ' 0 0 1 2 3 \n 0 1 2 3 4 \n 1 2 3 4 5 \n 1 2 3 4 0 \n 1 2 3 0 0'
    b = ' 1 \n 2 \n 3 \n 4 \n 5'
    x = ' -4 \n -1 \n 6 \n -2 \n -1'
    f = open('DBG/banded_A.txt', 'w')
    f.write(a)
    f.close();
    f = open('DBG/b.txt', 'w')
    f.write(b)
    f.close()
    f = open('DBG/x.txt', 'w')
    f.write(x)
    f.close();
    verify("DBG/banded_A.txt", "DBG/b.txt", "DBG/x.txt", l=2, u=2);

if(__name__=="__main__"):
    inlist = sys.argv
    if(len(inlist)==1): #tridiag
        verify("DBG/triA.txt", "DBG/trib.txt", "DBG/trix.txt");
    else: #banded 
        l,u = inlist[1], inlist[2]
        verify("DBG/banded_A.txt", "DBG/b.txt", "DBG/x.txt", \
                   l=int(l), u=int(u));
    #test1();



