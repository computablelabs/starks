import random
from starks.euclidean import gcd
from starks.euclidean import extended_euclidean_algorithm
from starks.polynomial import polynomials_over
from starks.poly_utils import generate_irreducible_polynomial
from starks.modp import IntegersModP
from starks.numbertype import FieldElement
from starks.numbertype import DomainElement
from starks.numbertype import memoize
from starks.numbertype import typecheck


@memoize
def FiniteField(p, m, polynomialModulus=None):
  """Create a type constructor for the finite field of order p^m for p prime, m >= 1"""
  Zp = IntegersModP(p)
  if m == 1:
    return Zp

  Polynomial = polynomials_over(Zp)
  if polynomialModulus is None:
    polynomialModulus = generate_irreducible_polynomial(modulus=p, degree=m)

  class Fq(FieldElement):
    field_size = int(p**m)
    primeSubfield = Zp
    ideal_generator = polynomialModulus
    operatorPrecedence = 3

    def __init__(self, poly):
      if isinstance(poly, bytes):
        self.poly = Polynomial(poly)
      elif type(poly) is Fq:
        self.poly = poly.poly
      elif type(poly) is int or type(poly) is Zp:
        self.poly = Polynomial([Zp(poly)])
      elif isinstance(poly, Polynomial):
        self.poly = poly % polynomialModulus
      else:
        self.poly = Polynomial([Zp(x) for x in poly]) % polynomialModulus

      self.field = Fq

    @typecheck
    def __add__(self, other):
      return Fq(self.poly + other.poly)

    @typecheck
    def __sub__(self, other):
      return Fq(self.poly - other.poly)

    @typecheck
    def __mul__(self, other):
      return Fq(self.poly * other.poly)

    @typecheck
    def __eq__(self, other):
      return isinstance(other, Fq) and self.poly == other.poly

    @typecheck
    def __ne__(self, other):
      return not self == other

    def __pow__(self, n):
      if n == 0: return Fq([1])
      if n == 1: return self
      if n % 2 == 0:
        sqrut = self**(n // 2)
        return sqrut * sqrut
      if n % 2 == 1: return (self**(n - 1)) * self

    #def __pow__(self, n): return Fq(pow(self.poly, n))
    def __neg__(self):
      return Fq(-self.poly)

    def __abs__(self):
      return abs(self.poly)

    def __repr__(self):
      # The code \u2208 prints \in in pretty symbols
      return repr(self.poly) + ' \u2208 ' + self.__class__.__name__

    @typecheck
    def __divmod__(self, divisor):
      q, r = divmod(self.poly, divisor.poly)
      return (Fq(q), Fq(r))

    def inverse(self):
      if self == Fq(0):
        raise ZeroDivisionError

      x, y, d = extended_euclidean_algorithm(self.poly, self.ideal_generator)
      if d.degree() != 0:
        raise Exception(
            'Somehow, this element has no inverse! Maybe intialized with a non-prime?'
        )

      return Fq(x) * Fq(d.coefficients[0].inverse())

    # TODO(rbharath): This function is broken!!
    def to_bytes(self):
      return self.poly.to_bytes()

  Fq.__name__ = 'F_{%d^%d}' % (p, m)
  Fq.p = p
  Fq.m = m
  Fq.base_field = Zp
  return Fq
