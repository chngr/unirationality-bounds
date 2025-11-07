"""
This file contains various functions that compute the numerical
functions appearing in the paper

  Raymond Cheng, "Unirationality of hypersurfaces via highly tangent lines".

Multi-degrees are represented here as multiplicity sequences, namely,
    mu[i] := number of times the degree i+1 appears.
"""
from math import comb, ceil
from fractions import Fraction

def degree_to_mu(d, m=1):
    """
    Return the multiplicity sequence
        mu = (0,...,0,m)
    corresponding to the multi-degree (d^m).
    """
    return [0 if i < d-1 else m for i in range(d)]

def is_linear(mu):
    """
    Check if a multiplicity sequence mu = (mu_1,...,mu_d) corresponds to
    a linear multi-degree, i.e. mu_i = 0 for all i > 1.
    """
    for i, m in enumerate(mu):
        if i > 0 and m > 0:
            return False
    return True

def is_quadric(mu):
    """
    Check if a multiplicity sequence mu = (mu_1,...,mu_d) corresponds to 
    a quadric hypersurface, i.e. mu_i = 0 for i > 2 and mu_2 = 1.
    """
    for i, m in enumerate(mu):
        if i == 1 and (m == 0 or m > 1):
            return False
        if i > 1 and m > 0:
            return False
    return True

def trim(mu):
    """
    Remove trailing zeroes in the multiplicity sequence mu.
    """
    while mu != [] and mu[-1] == 0:
        mu.pop()

def penta(mu):
    """
    Given a multiplicity sequence mu = (mu_1,...,mu_d), return the sequence
      mu' = (mu_1 + ... + mu_d,
             ...,
             mu_{d-2} + mu_{d-1} + mu_d,
             mu_{d-1} + mu_d - 1,
             mu_d - 1)
    corresponding to the penultimate tangent transform applied to mu. Note that
    mu' = 0 for the base case where mu corresponds to a linear multi-degree.
    """
    if is_linear(mu):
        return []
    else:
        trim(mu)
        mu_prime = [sum(mu[i:]) for i in range(len(mu))]
        mu_prime[-1] -= 1
        if len(mu) > 1:
            mu_prime[-2] -= 1
        return mu_prime

def pointed_lines(mu):
    """
    Given a multiplicity sequence mu = (mu_1,...,mu_d), return the sequence
        mu1 = (mu_1 + ... + mu_d,
               ...,
               mu_{d-2} + mu_{d-1} + mu_d,
               mu_{d-1} + mu_d - 1,
               mu_d - 1)
    corresponding to the pointed line construction applied to mu.
    """
    return [sum(mu[i:]) for i in range(len(mu))]


def r0(mu):
    """
    Given a multiplicity sequence mu, compute the quantity r0 from 1.7.
    """
    return sum([m * mu[m] for m in range(len(mu))]) - 1

def r(mu):
    """
    Given a multiplicity sequence mu, compute the function r from 1.12.
    """
    trim(mu)
    if mu == []:
      return -2
    r_max = r0(mu)
    mu_prime = mu
    depth = 0
    while not is_linear(mu_prime):
        mu_prime = penta(mu_prime)
        depth = depth + 1
        r_max = max(r_max, r0(mu_prime) + depth)
    return r_max

def n0(mu,r):
    """
    Given a multiplicity sequence mu and integer r, compute the ceiling
    of the quantity n0(mu,r) from p.9. The initial cases take care of those
    not covered in 1.11, and are the base cases defined in 1.12.
    """
    if r < -1:
        return 0
    elif r == -1:
        return sum(mu) - 1
    elif r == 0:
        return sum(pointed_lines(mu))
    elif is_quadric(mu):
        return 2*r + sum(mu) + 1
    else:
        return r + ceil(Fraction(1,r)*(sum([mu[d] * comb(d+1+r,r) for d in range(len(mu))]) - 1))

def n(mu,r):
    """
    Given a multiplicity sequence mu and integer r, compute the function
    n(mu,r) found in 1.12.
    """
    n_max = n0(mu,r)
    mu_prime = mu
    depth = 0
    while not is_linear(mu_prime):
        mu_prime = penta(mu_prime)
        depth = depth + 1
        n_max = max(n_max, n0(mu_prime,r-1) + 1)
    return n_max

def n_min(mu):
    """
    Given a multiplicity sequence mu, compute the function n(mu,r(mu)) given in
    1.14 for when a general complete intersection with multiplicity sequence mu
    becomes unirational.
    """
    return n(mu, r(mu))
