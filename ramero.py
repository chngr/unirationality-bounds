"""
This file contains numerical functions which compute the bounds contained in
Ramero's paper, "Effective estimates for unirationality".
"""
from math import ceil, comb
from fractions import Fraction

def e(ds, q):
    """
    Given a multi-degree ds = (d1,...,dc) and an integer q, compute
    the function e defined at the end of §3.
    """
    N = sum([comb(d + q, q) for d in ds])
    return q + ceil(Fraction(N,q+1))

def w(ds):
    """
    Given a multi-degree ds = (d1,...,dc), compute the function w
    defined at the beginning of §4.
    """
    if ds == []:
        return 0
    elif ds[0] == 1:
        return w(ds[1:])
    elif len(ds) == 1 and ds[0] == 2:
        return 0
    else:
        ds_prime = [d - 1 for d in ds]
        return max(r(ds_prime) - 1, e(ds_prime, w(ds_prime)))

def r(ds):
    """
    Given a multi-degree ds = (d1,...,dc), compute the function r
    defined at the beginning of §4.
    """
    if ds == []:
        return 0
    elif ds[0] == 1:
        return r(ds[1:]) + 1
    elif len(ds) == 1 and ds[0] == 2:
        return 2
    else:
        return 1 + sum([comb(d - 1 + w(ds), d - 1) - 1 for d in ds])

def m(d):
    """
    Given an integer d, compute the quantity m(d) appearing in Theorem 2.
    """
    ds = [d]
    return max(e(ds,w(ds)), r(ds))
