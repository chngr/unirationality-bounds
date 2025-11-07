import math
from fractions import Fraction

from sympy import binomial, fps, Mul, Pow, Rational, Integer, simplify, apart, together
from sympy.core.numbers import NegativeOne
from sympy.abc import x

class CoeffComputer:
    """
    An object class used to memoize computations of the quantities
      - F_i(x) as defined on p.11,
      - m_{i,j} as defined on p.11,
      - b_{i,j} as defined on p.13, and
      - c_{i,j} as defined on p.13.
    """

    def __init__(self):
        self._ms = {}
        self._cs = {}
        self._Fs = {0: 1}

    def m(self,i,j):
        """
        Compute the coefficients m_{i,j} of the power series
          F_i(x) = sum_{j ≥ 0} m_{i,j} x^{i+j}
        as defined on p.11. The computation is performed using the
        recursive formulae given in Lemma 2.6:

          m_{i+1,0} = (1/2) * m_{i,0}^2 - (1/2) * m_{i,0} + m_{i,1}, and
          m_{i+1,j} = 1/(j+2) * binomial(m_{i,0}+j-1,j) * (m_{i,0}^2 + (j-1)*m_{i,0} + 2)
                      + sum_{k = 0}^j binomial(m_{i,0}+j-k-1,j-k) * m_{i,k+1},
        
        where the base cases are m_{0,0} = 1 and m_{0,j} = 0 for all j ≥ 1.
        """
        if i == 0:
            return 1 if j == 0 else 0
        elif (i,j) in self._ms:
            return int(self._ms[(i,j)])
        else:
            m_prime = self.m(i-1,0)
            if j == 0:
                m = Rational(1,2) * (m_prime**2 - m_prime) + self.m(i-1,1)
                self._ms[(i,j)] = m
            else:
                m = Rational(1,j+2) * binomial(m_prime + j - 1, j) * \
                    (m_prime**2 + (j-1) * m_prime + 2) + \
                    sum([binomial(m_prime + j - k - 1, j - k) * self.m(i-1,k+1) for k in range(j+1)])
                self._ms[(i,j)] = m
            return int(m)

    def b(self, i, j):
        """
        Compute the quantity
          b_{i,j} := binomial(m_{i,0}+j-1,j} * m_{i,0}^j
        """
        return math.comb(self.m(i,0) + j - 1, j)/self.m(i,0)**j

    def c(self, i, j):
        """
        Compute the coefficients that ought to satisfy the bound
            m(i,j) <= c(i,j) * m(i,0) ** (1 + j/2)
        for each j >= 1 and when i >= 7.
        """
        if j == 0:
            raise Exception
        elif i <= 6:
            raise Exception
        elif i == 7:
            return 1
            return self.m(i,j)/self.m(i,0)**(1+j/2)
        elif (i,j) in self._cs:
            return self._cs[(i,j)]
        else: 
            r = self.m(i,0)
            c = 2**(1+j/2)*((self.b(i-1,j)/(j+2))*(1 + (j-1)/math.sqrt(2*r) + 1/r) \
                + sum([self.b(i-1,j-k)*self.c(i-1,k+1)/(2*r)**((k+1)/4) for k in range(j+1)]))
            self._cs[(i,j)] = c
            return c

    def F(self,i):
        if i not in self._Fs:
            self._Fs[i] = apart(Delta(i-1,self.F(i-1), m=self.m(i-1,0)) * x**(-i)) * x**i
        return self._Fs[i]

    def n(self,d):
        r = sum([self.m(i,0) for i in range(d-1)])
        return r + math.ceil(Fraction(math.comb(d+r,d) - 1,r))


def Delta(i,F,m=1):
    """
    Return the operator Delta_i^m(F) applied to a power series F.
    """
    return (1-x)**(-r) * F + ((1-x)**(-m) - 1) * (x**(i+1) - x**(i-1))

def coeff_formula(F):
    """
    Given an expression where F + 1 is a sum of terms of the form
        c * (x-1) **(-m)
    return the formula that extracts the k-th coefficient of the associated
    formal power series. We make the assumption that an expression that is
    of type
        - Pow looks like (x-1)**(-m)
        - Mult looks like c*(x-1)**(-m)
    """
    formula = 0
    for arg in apart(F).args:
        if arg.func == Pow:
            for arg2 in arg.args:
                if arg2.func == Integer or arg2.func == NegativeOne:
                    formula += binomial(x-arg2-1,-arg2-1)
        elif arg.func == Mul:
            c = 0
            m = 0
            for arg2 in arg.args:
                if arg2.func == Integer or arg2.func == NegativeOne:
                    c = arg2
                elif arg2.func == Pow:
                    for arg3 in arg2.args:
                        if arg3.func == Integer or arg3.func == NegativeOne:
                            m = arg3
            formula += (-1)**m * c * binomial(x-m-1,-m-1)
    return simplify(formula)

