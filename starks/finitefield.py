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


#this function generates a binary list of 2^index
def num_to_binary_list_pow(index):
  output_list = []

  for i in range(index):
    output_list.append(0)
  output_list.append(1)

  return output_list 

#binary addition of two binary strings from geeksforgeeks
def add_binary_nums(x, y): 
  max_len = max(len(x), len(y)) 
  
  x = x.zfill(max_len) 
  y = y.zfill(max_len)  
          
  # initialize the carry 
  carry = 0
  
  output_list = []

  # Traverse the string 
  for i in range(max_len - 1, -1, -1): 
    r = carry 
    r += 1 if x[i] == '1' else 0
    r += 1 if y[i] == '1' else 0
    output_list.append(1) if r % 2 == 1 else output_list.append(0)
    carry = 0 if r < 2 else 1     # Compute the carry. 
          
  if carry !=0 : output_list.append(1)
  
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

    # TODO(rbharath): This function is broken!!
    def to_bytes(self):
      return self.poly.to_bytes()

    #this function computes 2^x where x is a number in binary finite field representation and output is in in binary finite field representation
    def two_pow(self):
      num = Fq([1])
    
      for i in range(self.poly.degree()+1):
        if str(self.poly.coefficients[i])[0] == '1':
          l = num_to_binary_list_pow(i+1)
          num = num * Fq(l)
      
      return num



    #less than and equality operations between numbers in binary finite field representation
    def LT(self, other):
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


    #this function computes regular binary addition
    def Regular_Binary_Addition(self, other):
      x = ''
      y = ''

      for i in range(self.poly.degree()+1):
        x += str(self.poly.coefficients[i])[0]

      for i in range(other.poly.degree()+1):
        y += str(other.poly.coefficients[i])[0]

      list = add_binary_nums(x, y)

      return Fq(list)


  Fq.__name__ = 'F_{%d^%d}' % (p, m)
  Fq.p = p
  Fq.m = m
  Fq.base_field = Zp
  return Fq
