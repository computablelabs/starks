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
from starks.finitefield import FiniteField


#this function generates a binary list represents index
def num_to_binary_list(index):
  output_list = []

  while index > 0:
    output_list.append(index%2)
    index = index // 2

  return output_list

@memoize
def FloatingPoint(field):
  class Fp:
    def __init__(self, v, p, z, s):
      #these are four ellements of floating point numbers
      self.v = v
      self.p = p
      self.z = z
      self.s = s
   

    @typecheck
    def __eq__(self, other):
      return self.z == other.z or (self.v == other.v and self.p == other.p and self.z == other.z and self.s == other.s)


    #this is addition operation of two floating point numbers in binary finite field representation
    @typecheck
    def __add__(self, other):
      polysOver = polynomials_over(IntegersModP(field.p))

      if self.z == 1:
        v_c = other.v
        p_c = other.p
        z_c = other.z
        s_c = other.s
        output = Fp(v_c, p_c, z_c, s_c)
        return output

      if other.z == 1:
        v_c = self.v
        p_c = self.p
        z_c = self.z
        s_c = self.s
        output = Fp(v_c, p_c, z_c, s_c)
        return output

      a = self.v.LT(other.v)
      b = self.p.LT(other.p)

      if self.p == other.p:
        if a == 1:
          v_min = self.v
          v_max = other.v
          s_c = other.s
        else:
          v_min = other.v
          v_max = self.v
          s_c = self.s
      else:
        if b == 1:
          v_min = self.v
          v_max = other.v
          s_c = other.s
        else:
          v_min = other.v
          v_max = self.v
          s_c = self.s

      if b == 1:
        p_min = self.p
        p_max = other.p
      else:
        p_min = other.p
        p_max = self.p

      Delta = p_max - p_min

      c = Delta.LT(field(polysOver(num_to_binary_list(field.m))))

      if c == 0:
        v_c = v_max
        p_c = p_max
      else:
        v_c = v_max * Delta.two_pow()
        v_c = v_c + v_min
        p_c = p_min

      if v_c == field(polysOver([0])):
        z_c = 1
      else:
        z_c = 0

      output = Fp(v_c, p_c, z_c, s_c)
      return output

    #this is subtraction operation of two floating point numbers in binary finite field representation and it uses addition operation
    @typecheck
    def __sub__(self, other):
      other.s = other.s ^ 1
      output = self + other
      return output

    #this is division operation of two floating point numbers in binary finite field representation
    @typecheck
    def __truediv__(self, other):
      polysOver = polynomials_over(IntegersModP(field.p))

      if other.z == 1:
        raise ZeroDivisionError
      if other.v == field(polysOver([0])):
        raise ZeroDivisionError

      z_c = self.z | other.z

      s_c = self.s ^ other.s

      if z_c == 1:
        v_c = field(polysOver([0]))
      else:
        v_c = self.v / other.v

      if z_c == 1:
        p_c = field(polysOver([0]))
      else:
        p_c = self.p - other.p
        p_c = p_c + num_to_binary_list(self.v.poly.degree() + other.v.inverse().poly.degree() - (self.v.m - 1))

      output = Fp(v_c, p_c, z_c, s_c)
      return output

    #this is multiplication operation of two floating point numbers in binary finite field representation
    @typecheck
    def __mul__(self, other):
      polysOver = polynomials_over(IntegersModP(field.p))

      z_c = self.z | other.z

      s_c = self.s ^ other.s

      if z_c == 1:
        v_c = field(polysOver([0]))
      else:
        v_c = self.v * other.v

      if z_c == 1:
        p_c = field(polysOver([0]))
      else:
        p_c = self.p + other.p
        p_c = p_c + num_to_binary_list(self.v.poly.degree() + self.v.poly.degree() - (self.v.m - 1))
 
      output = Fp(v_c, p_c, z_c, s_c)
      return output


    #this function converts a floating point to fixed point
    def Float_to_Fix(self):
      polysOver = polynomials_over(IntegersModP(field.p))
      if self.z == 1:
        return field(polysOver([0])), 0

      result_of_power = self.p.two_pow()
      output = self.v * result_of_power
      return output, self.s


    #this function computes exponentiation based on regular binary addition
    def __pow__(self, other):
      fixed_of_Exp, sign = other.Float_to_Fix()

      polysOver = polynomials_over(IntegersModP(field.p))
      counter = field(polysOver([0]))
      output = Fp(field(polysOver([1])), field(polysOver([0])), 0, 0)
      while counter != fixed_of_Exp:
        output = output * self
        counter = counter.Regular_Binary_Addition(field(polysOver([1])))

      if other.s == 1:
        output = Fp(field(polysOver([1])), field(polysOver([0])), 0, 0)/output

      return output

  return Fp
