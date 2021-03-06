#! /usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function

import sympy
from sympy.abc import n, x

def intg(k, scale = 1):
    """
    return integral of T_{n+k} as a dictionary and as a function of n
    """
    d = dict()
    d[k+1] = 1/(2*(n+k+1))*scale
    d[k-1] = 1/(-2*(n+k-1))*scale
    return d
    
def iterate_once(d):
    klist = list(d.keys())
    ans = dict()
    for k in range(min(klist)-1,max(klist)+2):
        ans[k] = 0
    
    for k in klist:
        dd = intg(k, scale = d[k])
        for kk in dd.keys():
            ans[kk] += dd[kk]
      
    for k in ans.keys():
        ans[k] = sympy.factor(ans[k])
    return ans

def integral(k):
    """
    return coeffs of T_{n+k} in  k-fold integral of T_n as a dictionary 
    and as functions of n
    """
    d = intg(0)
    for i in range(0,k-1):
        d = iterate_once(d)
    return d

def coeff(k):
    """
    return the coeffs of T_n in the k-fold integral of 
    T_{n-k} ... T_{n+k} as functions of n
    """
    d = integral(k)
    ans = dict()
    for k in d.keys():
        ans[-k] = sympy.factor(d[k].subs(n, n-k))
    return ans

def show_dict(d):
    for k in sorted(d.keys()):
        if d[k] == 0:
            continue
        print('{0}\t--->\t{1}'.format(k, d[k]))

def show_derv1():
    u = 16*x**5 - 20*x**3 + 5*x;
    from sympy.abc import a, b
    from sympy import diff
    print(u)
    print(diff(u, x, 4) - (a*a+b*b)*diff(u, x, 2) + a*a*b*b*u)

def show_derv2():
    from sympy.abc import t
    k = sympy.symbols('k', integer=True)
    def diff(f):
        g = sympy.diff(f, t)/-sympy.sin(t)
        g = sympy.simplify(g)
        return g
    f = sympy.cos(k*t)
    g = diff(diff(f))
    print('------------------')
    print('2nd derv of T_k(x):\t{0}'.format(g))
    print('            at 0/1:\t{0}'.
          format(sympy.simplify(sympy.limit(g, t, 0))))
    print('          at pi/-1:\t{0}'.format(
        sympy.simplify(sympy.limit(g, t, sympy.pi))))
    g = diff(diff(g))
    print('4th derv of T_k(x):\t{0}'.format(g))
    print('            at 0/1:\t{0}'.
          format(sympy.simplify(sympy.limit(g, t, 0))))
    print('          at pi/-1:\t{0}'.format(
        sympy.simplify(sympy.limit(g, t, sympy.pi))))

if __name__ == '__main__':
    d = coeff(2)
    show_dict(d)
    print('----------------')
    d = coeff(4)
    show_dict(d)
    print('----------------')
    show_derv1()
    print('----------------')
    show_derv2()
    
