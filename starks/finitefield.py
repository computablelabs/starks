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

def num_to_binary_list_pow(index):
  output_list = []

  for i in range(index):
    output_list.append(0)
  output_list.append(1)

  return output_list 


def num_to_binary_list(index):
  output_list = []

  while index > 0:
    output_list.append(index%2)
    index = index // 2  

  return output_list


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
    
    # less than or equal operation
    def LQ(self, other):
      if self.poly.degree() > other.poly.degree():
        return 0
      elif other.poly.degree() > self.poly.degree():
        return 1

      xor_o = self + other

      if xor_o == Fq(0):
        return 1
      
      
      if str(self.poly.coefficients[xor_o.poly.degree()])[0] > str(other.poly.coefficients[xor_o.poly.degree()])[0]:
        return 0
      else:
        return 1

    # two pow self operation
    def pow(self):
      num = Fq([1])
    
      for i in range(self.poly.degree()+1):
        if str(self.poly.coefficients[i])[0] == '1':
          l = num_to_binary_list_pow(i+1)
          num = num * Fq(l)
      
      return num

    # addition floating point operation
    def addition(self, v_a, p_a, z_a, s_a, v_b, p_b, z_b, s_b):
      if z_a == 0:
        v_c = v_b
        p_c = p_b
        z_c = z_b
        s_c = s_b
        return v_c, p_c, z_c, s_c

      if z_b == 0:
        v_c = v_a
        p_c = p_a
        z_c = z_a
        s_c = s_a
        return v_c, p_c, z_c, s_c

      a = LQ(v_a.poly, v_b.poly)
      b = LQ(p_a.poly, p_b.poly)

      if p_a == p_b:
        if a == 1:
          v_min = v_a
          v_max = v_b
          s_c = s_b
        else:
          v_min = v_b
          v_max = v_a
          s_c = s_a
      else:
        if b == 1:
          v_min = v_a
          v_max = v_b
          s_c = s_b
        else:
          v_min = v_b
          v_max = v_a
          s_c = s_a

      Delta = p_max - p_min

      c = LQ(Delta.poly, num_to_binary_list(Fq.m))

      if c == 0:
        v_c = v_max
        p_c = p_max
      else:
        v_c = v_max * Delta.pow()
        v_c = v_c + v_min
        p_c = p_min

      if v_c == Fq(0):
        z_c = 1
      else:
        z_c = 0

      return v_c, p_c, z_c, s_c

    # division floating point operatin
    def division(self, v_a, p_a, z_a, s_a, v_b, p_b, z_b, s_b):
      z_c = z_a | z_b

      s_c = s_a ^ s_b

      if z_c == 1:
        v_c = Fq(0)
      else:
        v_c = v_a / v_b

      if z_c == 1:
        p_c = Fq(0)
      else:
        p_c = p_a - p_b
        p_c = p_c + num_to_binary_list(p_a.poly.degree() + p_b.inverse().poly.degree() - (Fq.m - 1))

      return v_c, p_c, z_c, s_c


    # multiplication floating point operation
    def multiplication(self, v_a, p_a, z_a, s_a, v_b, p_b, z_b, s_b):
      z_c = z_a | z_b

      s_c = s_a ^ s_b

      if z_c == 1:
        v_c = Fq(0)
      else:
        v_c = v_a * v_b

      if z_c == 1:
        p_c = Fq(0)
      else:
        p_c = p_a + p_b
        p_c = p_c + num_to_binary_list(p_a.poly.degree() + p_b.poly.degree() - (Fq.m - 1))

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
