try:
  from itertools import zip_longest
except ImportError:
  from itertools import izip_longest as zip_longest
import fractions

from starks.numbertype import memoize
from starks.numbertype import typecheck
from starks.numbertype import DomainElement
from starks.numbertype import Poly


def strip(L, elt):
  """strip all copies of elt from the end of the list"""
  if len(L) == 0: return L

  i = len(L) - 1
  while i >= 0 and L[i] == elt:
    i -= 1

  return L[:i + 1]

@memoize
def polynomials_over(ring=fractions.Fraction):
  """Create a polynomial with coefficients in a ring 
  
  Coefficients are in increasing order of monomial degree so that, for example,
  [1,2,3] corresponds to 1 + 2x + 3x^2.
  """

  class Polynomial(Poly):
    operatorPrecedence = 2

    @classmethod
    def factory(cls, L):
      return Polynomial([cls.ring(x) for x in L])

    def __init__(self, c):
      if isinstance(c, bytes):
        if len(c) % 32 != 0:
          raise ValueError("Bytelength must be multiple of 32")
        n_bytes = len(c) // 32
        coefficients = []
        for i in range(n_bytes):
          raw = c[32*i:32*(i+1)]
          coefficients.append(ring(raw))
        self.coefficients = coefficients

      elif type(c) is Polynomial:
        self.coefficients = c.coefficients
      elif isinstance(c, ring):
        self.coefficients = [c]
      elif not hasattr(c, '__iter__') and not hasattr(c, 'iter'):
        self.coefficients = [ring(c)]
      else:
        self.coefficients = [ring(coeff) for coeff in c]

      self.coefficients = strip(self.coefficients, ring(0))

    def is_zero(self):
      return self.coefficients == []

    def __repr__(self):
      if self.is_zero():
        return '0'

      return ' + '.join([
          '%s x^%d' % (a, i) if i > 0 else '%s' % a
          for i, a in enumerate(self.coefficients)
      ])

    def __abs__(self):
      return len(
          self.coefficients
      )  # the valuation only gives 0 to the zero polynomial, i.e. 1+degree

    def __len__(self):
      return len(self.coefficients)

    def __sub__(self, other):
      return self + (-other)

    def __iter__(self):
      return iter(self.coefficients)

    def __neg__(self):
      return Polynomial([-a for a in self])

    def iter(self):
      return self.__iter__()

    def leading_coefficient(self):
      return self.coefficients[-1]

    def degree(self):
      return abs(self) - 1

    @typecheck
    def __eq__(self, other):
      return self.degree() == other.degree() and all(
          [x == y for (x, y) in zip(self, other)])

    @typecheck
    def __ne__(self, other):
      return self.degree() != other.degree() or any(
          [x != y for (x, y) in zip(self, other)])

    @typecheck
    def __add__(self, other):
      new_coefficients = [
          sum(x) for x in zip_longest(self, other, fillvalue=self.ring(0))
      ]
      return Polynomial(new_coefficients)

    @typecheck
    def __mul__(self, other):
      if self.is_zero() or other.is_zero():
        return Zero()

      newCoeffs = [self.ring(0) for _ in range(len(self) + len(other) - 1)]

      for i, a in enumerate(self):
        for j, b in enumerate(other):
          newCoeffs[i + j] += a * b

      return Polynomial(newCoeffs)

    @typecheck
    def __divmod__(self, divisor):
      quotient, remainder = Zero(), self
      divisorDeg = divisor.degree()
      divisorLC = divisor.leading_coefficient()

      while remainder.degree() >= divisorDeg:
        monomialExponent = remainder.degree() - divisorDeg
        monomialZeros = [self.ring(0) for _ in range(monomialExponent)]
        monomialDivisor = Polynomial(
            monomialZeros + [remainder.leading_coefficient() / divisorLC])

        quotient += monomialDivisor
        remainder -= monomialDivisor * divisor

      return quotient, remainder

    @typecheck
    def __truediv__(self, divisor):
      if divisor.is_zero():
        raise ZeroDivisionError
      return divmod(self, divisor)[0]

    @typecheck
    def __mod__(self, divisor):
      if divisor.is_zero():
        raise ZeroDivisionError
      return divmod(self, divisor)[1]

    # TODO(rbharath): Possibly type-check this.
    def __call__(self, x):
      y = self.ring(0)
      power_of_x = 1
      for i, a in enumerate(self):
        y += power_of_x * a 
        power_of_x = (power_of_x * x)
      return y

    def to_bytes(self):
      coeff_bytes = [coeff.to_bytes(32, 'big') for coeff in self.coefficients]


  def Zero():
    return Polynomial([])

  Polynomial.ring = ring 
  Polynomial.__name__ = '(%s)[x]' % ring.__name__
  Polynomial.englishName = 'Polynomials in one variable over %s' % ring.__name__
  return Polynomial
