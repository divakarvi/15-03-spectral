#! /usr/bin/env python
import matplotlib as mpl
mpl.rcParams['font.size'] = 18
from matplotlib import pyplot as plt
from numpy import genfromtxt
fig = plt.figure(1)
ax = plt.subplot(111)

x = genfromtxt('sinpiya2e4_x0.dat', dtype='float')
y = genfromtxt('sinpiya2e4_y0.dat', dtype='float')
l, = plt.loglog(x,y,'k')
l.set_linewidth(2)
h1 = l

x = genfromtxt('sinpiya2e4_x1.dat', dtype='float')
y = genfromtxt('sinpiya2e4_y1.dat', dtype='float')
l, = plt.loglog(x,y,'k')
l.set_linestyle('--')
l.set_linewidth(2)
h2 = l

x = genfromtxt('sinpiya2e4_x2.dat', dtype='float')
y = genfromtxt('sinpiya2e4_y2.dat', dtype='float')
l, = plt.loglog(x,y,'k')
l.set_linestyle('-.')
l.set_linewidth(2)
h3 = l

plt.legend([h1, h2, h3], ['u error', 'du error', 'num diff error'], loc=2)
plt.xlabel('M')

ax.axis([3.0000e+01, 2.0000e+04, 1.0000e-16, 1.0000e-04]) 
#ax.set_title('a = 2e4, sin(pi*y)/single grid', fontsize=20) 
#plt.show() 
plt.savefig('sinpiya2e4.pdf')
