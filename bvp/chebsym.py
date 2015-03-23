from __future__ import absolute_import
from __future__ import print_function

import sympy
from sympy.abc import n, x

def intg(k, scale = 1):
    """
    return integral of T_{n+k} as a dictionary
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
    return k-fold integral of T_n as a dictionary
    """
    d = intg(0)
    for i in range(0,k-1):
        d = iterate_once(d)
    return d

def coeff(k):
    """
    return the coeff of T_n in the k-fold integral of a cheb series
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



if __name__ == '__main__':
    d = integral(4)
    show_dict(d)
    print('----------------')
    d = coeff(4)
    show_dict(d)
    
    
    
