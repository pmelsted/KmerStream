#!/usr/bin/env python
import sys
from math import exp
from scipy.optimize import brentq

def EMfit2(F0, f1, F1, k):
    epstol = 1e-8
    x = float(F1) / F0
    e = 0
    e0 = 1
    x0 = 0
    maxit = 1000
    it = 0
    while abs(x - x0) / x + abs(e - e0) > epstol and it < maxit:
        it = it + 1
        #print x, e, abs(x-x0)/x
        x0 = x
        e0 = e
        x1 = -100
        while abs(x - x1) / x > epstol / 2:
            x1 = x
            x  = float(F1) / F0 * (3 * k * (1 - exp(-x * e / (3 * k))) + 1 - exp(-x * (1 - e)))
        def func(t):
            return float(f1) / F1 - (t * exp(-x * t / (3 * k)) + (1 - t) * exp(-x * (1 - t)))

        e = brentq(func, 0, 1)
    return x, e

if __name__ == '__main__':
    v = sys.argv
    if len(v) <= 1:
        print("Usage: estimate.py TSV-FILE")
        exit(1)
    else:
        fn = v[1]
        with open(fn) as f:
            header = f.readline().strip()
            print(header + "\tF0-f1\tG\te\tlambda")
            for line in f:
                l = line.split()
                q,k = l[:2]
                F0,f1,F1 = [float(x) for x in l[2:]]
                x,e = EMfit2(F0,f1,F1,int(k))
                print(line.strip() +'\t'+ '\t'.join(str(t) for t in (F0-f1,F1/x,e,x)))
