"""Classes and functions to aid representatin of multivariate polynomials.

Multivariate polynomials are used to represent transitions between adjacent
computation states.
"""
from __future__ import annotations
from typing import List 
from typing import Tuple
from typing import Dict
from typing import Callable
from typing import Any
#from starks.poly_utils import construct_multivariate_coefficients
from starks.numbertype import memoize
from starks.numbertype import Field
from starks.numbertype import FieldElement
from starks.numbertype import MultiVarPoly
from starks.numbertype import typecheck 
from sympy import Poly
from sympy import div
from sympy import invert
from starks.modp import IntegersModP


def remove_zero_coefficients(coefficients: Dict) -> Dict:
  reduced = {}
  for monomial in coefficients:
    if coefficients[monomial] == 0:
      continue
    else:
      reduced[monomial] = coefficients[monomial]
  return reduced

def add_power_tuples(a: Tuple, b: Tuple) -> Tuple:
  add_l = []
  min_a_b = min(len(a), len(b))
  for i in range(min_a_b):
    add_l.append(a[i]+b[i])
   
  if len(a) == min_a_b:
    for i in range(len(b)-min_a_b):
      add_l.append(b[i])
  else:
    for i in range(len(a)-min_a_b):
      add_l.append(a[i])
  return tuple(add_l)

def sum_power_tuple(a: Tuple) -> Any:
  out = 0
  for a_i in a:
    out += a_i
  return out

# TODO(rbharath): How does the memoization code actually work?
@memoize
def multivariates_over(ring: Field, num_vars: int) -> MultiVarPoly:
  """Create a multivariate polynomial.

  Let R be the ring and n = num_vars. Then the polynomial ring we're
  constructing here is R[X_1,...,X_n].
  """

  class MultivariatePolynomial(MultiVarPoly):
    # TODO(rbharath): This operator precedence bit is black magic. This needs
    # to be handled more systematically.
    #operatorPrecedence = 2
    operatorPrecedence = 4

    # TODO(rbharath): Using Any here isn't optimal. cls is meant to be a ring type. 
    @classmethod
    def factory(cls: Any, coefficients: Dict = None,
                step_fn: Callable = None) -> MultivariatePolynomial:
      """Constructs a multivariate polynomial with given coefficients."""
      if coefficients is not None:
        return MultivariatePolynomial(coefficients) 
      elif step_fn is not None:
        coefficients = construct_multivariate_coefficients(step_fn)
        return MultivariatePolynomial(coefficients)

    def __init__(self, c):
      if type(c) is MultivariatePolynomial:
        self.coefficients = c.coefficients
      elif isinstance(c, ring):
        self.coefficients = {(0,)*num_vars: c}
      elif isinstance(c, dict):
        coefficients = remove_zero_coefficients(c)
        self.coefficients = c
      elif isinstance(c, int):
        self.coefficients = {(0,)*num_vars: ring(c)}
      else:
        raise ValueError
      self.coefficients = remove_zero_coefficients(self.coefficients)

    def __len__(self):
      return len(self.coefficients)

    def is_zero(self):
      return self.coefficients == {}

    def __repr__(self):
      if self.is_zero():
        return '0'

      def power_tuple_to_string(power_tup):
        return "".join(["*X_%d**%d" % (i+1, power) for (i, power) in enumerate(power_tup)])

      return ' + '.join([
          '%s %s' % (str(coeff), power_tuple_to_string(power_tup)) if power_tup != (0,)*num_vars else '%s' % coeff 
          for power_tup, coeff in self.coefficients.items()
      ])

    def degree(self):
      """TODO(rbharath): Computing the degree is a little tricky."""
      max_deg = 0
      for power_tup in self.coefficients.keys():
        if sum_power_tuple(power_tup) > max_deg:
          max_deg = sum_power_tuple(power_tup)
      return max_deg

    def __sub__(self, other):
      return self + (-other)

    def __iter__(self):
      # Tuples are sorted in dictionary orderr
      keys = sorted(self.coefficients.keys())
      for key in keys:
        yield (key, self.coefficients[key])

    def __neg__(self):
      return MultivariatePolynomial({power_tup: -coeff for (power_tup, coeff) in self})

    @typecheck
    def __eq__(self, other):
      return (self.degree() == other.degree() and 
          # Checks monomials are the same
          self.coefficients.keys() == other.coefficients.keys() and
          all(
          [self[key] == other[key] for key in self.coefficients.keys()]))

    def __getitem__(self, power_tup: Tuple[int, ...]) -> FieldElement:
      # Monomials not present in multivariate polynomial have coefficient 0
      if power_tup in self.coefficients:
        return self.coefficients[power_tup]
      else:
        return ring(0)

    @typecheck
    def __add__(self, other):
      self_monomials = set(self.coefficients.keys())
      other_monomials = set(other.coefficients.keys())
      joint_monomials = self_monomials.union(other_monomials)
      new_coefficients = { 
          monomial: self[monomial] + other[monomial] for monomial in joint_monomials}
      return MultivariatePolynomial(new_coefficients)

    @typecheck
    def __mul__(self, other):
      if self.is_zero() or other.is_zero():
        return Zero()

      new_coeffs = {}
      for i, (a, a_coeff) in enumerate(self):
        for j, (b, b_coeff) in enumerate(other):
          prod = add_power_tuples(a, b)
          coeff = a_coeff * b_coeff 
          if prod not in new_coeffs:
            new_coeffs[prod] = ring(0)
          new_coeffs[prod] += coeff 
      return MultivariatePolynomial(new_coeffs)

    @typecheck
    def __truediv__(self, divisor):
      """" 
      div operation is implemented by using div in sympy
      The approach is, we convert MultivariatePolynomial to string which will be used to build MultivariatePolynomial in sympy,
      then we use division which is already available in sympy, the output of division then converted to MultivariatePolynomial
      which is compatible with our design
      """

      X = Poly(str(self))
      X_size = self.size_p()
      Y = Poly(str(divisor))
      Y_size = divisor.size_p()
      Z = div(X, Y)

      # the output of the division converted again to MultivariatePolynomial
      Z_str = str(Z)
      result = ""
      i = 6
      while i < len(Z_str) and Z_str[i] != ",":
        result += Z_str[i]
        i += 1

      # now we have string format of the output, it should be converted to MultivariatePolynomial
      return self.StrToMulti(result, max(X_size, Y_size))

    # it returns maximum degree in a polynomial
    def size_p(self):
      Max_Sym = 0
      st = str(self)
      for i in range(len(st)):
        if st[i] == "X" and i+2 < len(st) and  st[i+1] == "_" and int(st[i+2]) > Max_Sym:
          Max_Sym = int(st[i+2])

      return Max_Sym
    
    # concvert the string to MultivariatePolynomial, the approach is that all the  should be extracted to define MultivariatePolynomial agaian 
    # it works as a parser
    def StrToMulti(self, st, s):
      Max_Sym = 0
      for i in range(len(st)):
        if st[i] == "_" and int(st[i+1]) > Max_Sym:
          Max_Sym = int(st[i+1])

      zero = []
      for i in range(max(s, Max_Sym)):
        zero.append(0)
      zero_t = tuple(zero)
      new_coeffs = {}

      if Max_Sym == 0:
        new_coeffs[zero_t] = ring(0)
        new_coeffs[zero_t] += int(st)
        return MultivariatePolynomial(new_coeffs)

      i = 0
      while i < len(st):
        temp = ""
        while st[i].isdigit():
          temp += st[i]
          i += 1

        if temp == "":
          temp_i = 1
        else:
          temp_i = int(temp)

        temp_zero = zero_t

        while i < len(st) and st[i] is not " ":
          if st[i] == "X" and st[i+1] == "_":
            i += 2
          elif st[i] == "*" and st[i+1] == "X" and st[i+2] == "_":
            i += 3

          index = ""
          while i < len(st) and st[i].isdigit(): 
            index += st[i]
            i += 1
          index_i = int(index)

          if i >= len(st):
            pow_i = 1
          else:
            pow_i = 1
            if st[i] == "*" and i+1 < len(st) and st[i+1] == "*":
              i += 2
              poww = ""
              while st[i].isdigit():
                poww += st[i]
                i += 1
              pow_i = int(poww)

            elif st[i] == "*" and i+1 < len(st) and st[i+1] is not "*":
              pow_i = 1

          temp_zero_l = list(temp_zero)
          temp_zero_l[index_i-1] = pow_i
          temp_zero = tuple(temp_zero_l)

        new_coeffs[temp_zero] = ring(0)
        new_coeffs[temp_zero] += temp_i
        i += 3

      return MultivariatePolynomial(new_coeffs) 

    # compute remainder of a division
    def div_remainder(self, Y, Z):
      return self-Y*Z   

    # we use div implementation in sympy, this implementation is different from truediv because we may have
    # two arguments for the div that have differrent types, one MultivariatePolynomial and one polynomial with one variable
    # so we need to check the type of both arguments and make sure both are MultivariatePolynomial before division
    def division(self, divisor):
      X = Poly(self.CheckforDiv(self))
      X_size = self.size_p()
      Y = Poly(self.CheckforDiv(divisor))
      Y_size = 1      

      Z = div(X, Y)
      print(Z)
      Z_str = str(Z)
      result = ""
      i = 6
      while i < len(Z_str) and Z_str[i] != ",":
        result += Z_str[i]
        i += 1

      return self.StrToMulti(result, max(X_size, Y_size))

    # both elements of division should be in the same type, this fuction converts a polynomial with one variable to MultivariatePolynomial
    def PolytoMulti(self, p):
      size = self.size_p()

      Y = self.CheckforDiv(p)

      return self.StrToMulti(Y, size)
 
    # check the elements in division has the correct type, if it is polynomial with one varibale, we need to make it ready to convert it to MultivariatePolynomial
    def CheckforDiv(self, p):
      poly = str(p)
      if poly.find("F_") == -1:
        return poly

      result = ""
      i = 0
      while i < len(poly):
        if i+4 < len(poly):
          result += poly[i]
          if poly[i+4] == "F":
            while poly[i] is not "}":
              i += 1
        else:
          result += poly[i]
        i += 1

      return result         

    # TODO(rbharath): Possibly type-check this.
    def __call__(self, vals):
      assert len(vals) == num_vars
      y = ring(0)
      power_of_x = 1
      for _, (a, a_coeff) in enumerate(self):
        prod = ring(1)
        for i, power in enumerate(a):
          prod *= vals[i]**power
        y += a_coeff * prod
      return y

  def Zero():
    return MultivariatePolynomial({})


  MultivariatePolynomial.ring = ring 
  MultivariatePolynomial.num_vars = num_vars
  MultivariatePolynomial.__name__ = "".join(["(%s)" % ring.__name__, "[", ",".join(["X_%d" % (i+1) for i in range(num_vars)]), "]"])
  return MultivariatePolynomial
