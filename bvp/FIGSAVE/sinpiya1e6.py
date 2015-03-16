#! /usr/bin/env python
import matplotlib as mpl
mpl.rcParams['font.size'] = 18
from matplotlib import pyplot as plt
from numpy import genfromtxt

fig = plt.figure(1)
ax = plt.subplot(111)

x = genfromtxt('sinpiya1e6_x0.dat', dtype='float')
y = genfromtxt('sinpiya1e6_y0.dat', dtype='float')
l, = plt.loglog(x,y, 'k')
l.set_linewidth(2)
h1 = l

x = genfromtxt('sinpiya1e6_x1.dat', dtype='float')
y = genfromtxt('sinpiya1e6_y1.dat', dtype='float')
l, = plt.loglog(x,y, 'k')
l.set_linestyle('--')
l.set_linewidth(2)
h2 = l

x = genfromtxt('sinpiya1e6_x2.dat', dtype='float')
y = genfromtxt('sinpiya1e6_y2.dat', dtype='float')
l, = plt.loglog(x,y, 'k')
l.set_linestyle('-.')
l.set_linewidth(2)
h3 = l

plt.legend([h1, h2, h3], ['u error', 'du error', 'num diff error'], loc=2)
plt.text(0.05, 0.5, '$a=10^6$', fontsize = 18, transform = plt.gca().transAxes)
ax.axis([3.0000e+01, 2.0000e+04, 1.0000e-16, 1.0000e-04]) 
#plt.show() 
plt.xlabel('M')
plt.savefig('sinpiya1e6.pdf')
