# This file verifies various numerical assertions made in the accompanying paper
import math
import sympy

import ramero
from dimensions import *
from power_series import CoeffComputer

# Initialize the object to compute the coefficients m(i,j)
C = CoeffComputer()

def verify_introduction():
    # Verify the two numerical assertions in the introduction
    assert(C.n(10) == 192884152577980851363553858004926940342106493833715693762179)
    assert(196 < math.log2(C.n(10)) < 197)
    assert(171550 < math.log2(ramero.m(10)) < 171551)

def verify_values():
    # Check the numerical values n(d) for 3 ≤ d ≤ 9 given in the beginning of §2
    values = [4, 9, 22, 160, 20376, 11914188890, 8616199237736295920955120]
    for d in range(3,10):
        assert(n_min(degree_to_mu(d)) == values[d-3])

    # Check the values given in Figure 1 on p.12
    table = [[1,        3,           4,               5],
             [3,        8,           13,              19],
             [11,       48,          127,             275],
             [103,      1106,        7051,            33955],
             [6359,     485280,      21029990,        654279500],
             [20700541, 88819638509, 214404499562520, 368104651084030885]]
    for i in range(6):
        for j in range(4):
            assert(C.m(i+3,j) == table[i][j])
    return True

def verify_2_7():
    # Check the inequality m(i,0)^2 < 2*m(i+1,0) for 1 ≤ i ≤ 4 for Lemma 2.7.
    for i in range(1,5):
        assert(C.m(i,0)**2 < 2*C.m(i+1,0))
    return True

def verify_2_9():
    # Check the inequality m(7,1) < m(7,0)^(3/2) and m(7,2) < m(7,0)^2 for the
    # beginning of the proof of 2.9.
    assert(C.m(7,1) < C.m(7,0)**(3/2))
    assert(C.m(7,2) < C.m(7,0)**2)

    # Verify the inequality between log(m(7,0))/2 and the harmonic-like number
    #     sum_{l = 0}^{119} 1/(2 + l)
    # found in the middle of 2.9.
    assert(math.log(C.m(7,0))/2 > sum([1/(2+l) for l in range(120)]))

    # Verify the explicit bounds on the b(7,j) for 0 ≤ j ≤ 4 in 2.9
    assert(C.b(7,0) == 1)
    assert(C.b(7,1) == 1)
    assert(C.b(7,2) < 2/3)
    assert(C.b(7,3) < 1/4)
    assert(C.b(7,4) < 1/16)

    # Verify the bound (*) in 2.9 for j = 1 and j = 2
    assert((1/3) * (1+1/C.m(8,0)) +
      sum([C.b(7,1-k) * (2*C.m(8,0))**(-(k+1)/4) for k in range(2)]) < 2**(-3/2))
    assert((C.b(7,2)/4) * (1 + (2*C.m(8,0))**(-1/2) + 1/C.m(8,0)) +
      sum([C.b(7,2-k) * (2*C.m(8,0))**(-(k+1)/4) for k in range(3)]) < 2**(-5/2))

    # Verify the explicit inequality at the end of 2.9.
    assert((1/5) * (1 + 2/4**6 + 1/4**12) + 1/60 + 2/(3 * 4**5) + 1/4**8 + 1/4**11 < 1/4)
    return True

def verify_2_12():
    # Verify 2.12 for i = 6
    assert(sum([C.m(k,0) for k in range(7)]) <= 2**(2**3))
    return True

def verify_2_13():
    # Verify the polynomial given in 2.13
    ms = sympy.symbols("m:5")
    x = sympy.Symbol("x")
    expr1 = sympy.expand_func(sum([ms[d] * (sympy.binomial(d+x,x) - x*d - 1) for d in range(1,5)])).simplify()
    expr2 = (x/24) * (-(12*ms[2] + 28*ms[3] + 46*ms[4])
                      +(12*ms[2] + 24*ms[3] + 35*ms[4]) * x
                      +(4*ms[3] + 10*ms[4]) * x**2
                      + ms[4] * x**3).expand()
    assert(expr1 == expr2)
    return True

def verify_2_14():
    # Verify the statement for 3 ≤ d ≤ 7
    for d in range(3,8):
        mu = degree_to_mu(d)
        R = sum([C.m(i,0) for i in range(d-1)])
        assert(n_min(mu) == n0(mu,R))

    # Verify the inequality m(d-4,3) ≤ m(d-4,0) + m(d-3,0) + m(d-2,0) - 2d - 1
    # in Step 1 for 8 ≤ d ≤ 10.
    for d in range(8,11):
        assert(C.m(d-4,3) <= C.m(d-4,0) + C.m(d-3,0) + C.m(d-2,0) - 2*d - 1)

    # Verify the explicit inequality for d = 11 in Step 1.
    d = 11
    assert(4*C.m(d-2,0)**(5/8) <= C.m(d-4,0) + C.m(d-3,0) + C.m(d-2,0) - 2*d - 1)
    return True

def verify_2_15():
    # Verify the statement for 6 ≤ d ≤ 7
    for d in range(6,8):
        assert(n_min(degree_to_mu(d)) <= 2**((d-1)*2**(d-5)))
    return True

verify_introduction()
verify_values()
verify_2_7()
verify_2_9()
verify_2_12()
verify_2_13()
verify_2_14()
verify_2_15()
