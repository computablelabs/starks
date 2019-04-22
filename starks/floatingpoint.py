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

#this function generates a binary list of 2^index
def num_to_binary_list_pow(index):
  output_list = []

  for i in range(index):
    output_list.append(0)
  output_list.append(1)

  return output_list 


#this function generates a binary list represents index
def num_to_binary_list(index):
  output_list = []

  while index > 0:
    output_list.append(index%2)
    index = index // 2  

  return output_list


#less than and equality operations between numbers in binary finite field representation
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


#this function computes 2^self where self is a number in binary finite field representation and output is in in binary finite field representation
def pow(self):
  num = Fq([1])
    
  for i in range(self.poly.degree()+1):
    if str(self.poly.coefficients[i])[0] == '1':
      l = num_to_binary_list_pow(i+1)
      num = num * Fq(l)
      
  return num


@memoize
def FloatingPoint(field):
  class Fp:
    def __init__(self, v, p, z, s):
      #these are four ellements of floating point numbers
      self.v = v
      self.p = p
      self.z = z
      self.s = s
   

    #this is addition operation of two floating point numbers in binary finite field representation
    @typecheck
    def __add__(self, other):
      if self.z == 0:
        v_c = other.v
        p_c = other.p
        z_c = other.z
        s_c = other.s
        output = Fp(v_c, p_c, z_c, s_c)
        return output

      if other.z == 0:
        v_c = self.v
        p_c = self.p
        z_c = self.z
        s_c = self.s
        output = Fp(v_c, p_c, z_c, s_c)
        return output

      a = LQ(v_a.poly, v_b.poly)
      b = LQ(p_a.poly, p_b.poly)

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

      output = FloatingPointNumber(v_c, p_c, z_c, s_c)
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
      z_c = self.z | other.z

      s_c = self.s ^ other.s

      if z_c == 1:
        v_c = Fq(0)
      else:
        v_c = self.v / other.v

      if z_c == 1:
        p_c = Fq(0)
      else:
        p_c = self.p - other.p
        p_c = p_c + num_to_binary_list(self.p.poly.degree() + other.p.inverse().poly.degree() - (self.v.m - 1))

      output = Fp(v_c, p_c, z_c, s_c)
      return output

    #this is multiplication operation of two floating point numbers in binary finite field representation
    @typecheck
    def __mul__(self, other):
      z_c = self.z | other.z

      s_c = self.s ^ other.s

      if z_c == 1:
        v_c = Fq(0)
      else:
        v_c = self.v * other.v

      if z_c == 1:
        p_c = Fq(0)
      else:
        p_c = self.p + other.p
        p_c = p_c + num_to_binary_list(self.p.poly.degree() + self.p.poly.degree() - (self.v.m - 1))
 
      output = Fp(v_c, p_c, z_c, s_c)
      return output

  return Fp
