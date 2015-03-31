import numpy as np
import sympy as sym
from sympy.abc import x, a, b

def d(k):
    return 1.0/np.sqrt(2.0*(2*k+3)**2*(2*k+5))

def e(k): 
    return 2.0/(2.0*k+1)

def g(k):
    return (2.0*k+3)/(2.0*k+7)

def h(k):
    return -(1.0+g(k))

def bkk(k):
    return d(k)**2*(e(k) + h(k)**2*e(k+2) + g(k)**2*e(k+4))

def bkp2(k):
    return d(k)*d(k+2)*(h(k)*e(k+2) + g(k)*h(k+2)*e(k+4))

def bkp4(k):
    return d(k)*d(k+4)*g(k)*e(k+4)

def ckk(k):
    return -2.0*(2.0*k+3)*d(k)**2*h(k)
    
def ckp2(k):
    return -2.0*(2.0*k+2)*d(k)*d(k+2)

def gg(k):
    return -2.0*(2.0*k+5)/(2.0*k+7)

def L(x, n):
    assert isinstance(x, sym.Symbol)
    if n == 0:
        return 1
    elif n == 1:
        return x
    elif n > 1 and isinstance(n, int):
        return sym.simplify(((2*n-1)*x*L(x, n-1) -(n-1)*L(x,n-2))/n)

def LL(x, k):
    ans = L(x, k) + L(x, k+2)*(-2)*(2*k+5)/(2*k+7) + L(x, k+4)*(2*k+3)/(2*k+7)
    #ans = ans/sym.sqrt(2*(2*k+3)**2*(2*k+5))
    #ans = sym.simplify(ans)
    return ans

u1 = LL(x, 0)
f1 = sym.diff(u1, x, 4) - (a**2 + b**2)*sym.diff(u1, x, 2) + a**2*b**2*u1
f1 = sym.simplify(f1)

u2 = LL(x, 1)
f2 = sym.diff(u2, x, 4) - (a**2 + b**2)*sym.diff(u2, x, 2) + a**2*b**2*u2
f2 = sym.simplify(f2)
