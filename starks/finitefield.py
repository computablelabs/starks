import random
from starks.euclidean import gcd
from starks.euclidean import extended_euclidean_algorithm
from starks.polynomial import polynomials_over
from starks.poly_utils import generate_irreducible_polynomial
from starks.poly_utils import generate_primitive_polynomial
from starks.modp import IntegersModP
from starks.numbertype import FieldElement
from starks.numbertype import DomainElement
from starks.numbertype import memoize
from starks.numbertype import typecheck
from starks.poly_utils import is_irreducible


@memoize
def FiniteField(p, m, polynomialModulus=None):
  """Create a type constructor for the finite field of order p^m for p prime, m >= 1"""
  if polynomialModulus is not None:
    if not is_irreducible(polynomialModulus, p):
      raise ValueError("Must provide an irreducible polynomial as modulus.")
  Zp = IntegersModP(p)
  if m == 1:
    return Zp

  Polynomial = polynomials_over(Zp)
  if polynomialModulus is None:
    polynomialModulus = generate_primitive_polynomial(modulus=p, degree=m)

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
    
    def num_to_poly(n):
      if x == 0: return zero()
      if x == 1: return Fq([1])
      if x == 2: return Fq([0,1])

      one = Fq([1])
      two = Fq([0,1])
   
      num = Fq(Zero())

      if n%2 == 1:
        num = num + one

      n = n // 2
      count = 1
      while n > 0:
        two = two ** count
        count = count + 1
        if n%2 == 1:
          num = num + two
        n = n // 2

      return num

    def division(v_a, p_a, z_a, s_a, v_b, p_b, z_b, s_b):
      z_c = z_a | z_b

      s_c = s_a ^ s_b

      if z_c == 1:
        v_c = Fq(Zero())
      else:
        v_c = v_a / v_b

      if z_c == 1:
        p_c = Fq(Zero())
      else:
        p_c = p_a - p_b - num_to_poly((p_a.poly).degree - (p_b.poly).degree - (Fq.m - 1))

      return v_c, p_c, z_c, s_c


    def multiplication(v_a, p_a, z_a, s_a, v_b, p_b, z_b, s_b):
      z_c = z_a | z_b

      s_c = s_a ^ s_b

      if z_c == 1:
        v_c = Fq(Zero())
      else:
        v_c = v_a * v_b

      if z_c == 1:
        p_c = Fq(Zero())
      else:
        p_c = p_a + p_b + num_to_poly((p_a.poly).degree + (p_b.poly).degree - (Fq.m - 1)) 

      return v_c, p_c, z_c, s_c

    # TODO(rbharath): This function is broken!!
    def to_bytes(self):
      return self.poly.to_bytes()

    # TODO(rbharath): This function is broken!!
    def to_bytes(self):
      return self.poly.to_bytes()

  Fq.__name__ = 'F_{%d^%d}' % (p, m)
  Fq.p = p
  Fq.m = m
  Fq.base_field = Zp
  return Fq
