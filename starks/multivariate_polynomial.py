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

def remove_zero_coefficients(coefficients: Dict) -> Dict:
  reduced = {}
  for monomial in coefficients:
    if coefficients[monomial] == 0:
      continue
    else:
      reduced[monomial] = coefficients[monomial]
  return reduced

def add_power_tuples(a: Tuple, b: Tuple) -> Tuple:
  if len(a) != len(b):
    raise ValueError("Can't add tuples of different lengths")
  return tuple([a_i + b_i for a_i, b_i in zip(a, b)])

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
        return "".join(["X_%d^%d" % (i+1, power) for (i, power) in enumerate(power_tup)])

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
    def __divmod__(self, divisor):
      """"TODO(rbharath): Implementing polynomial division in multiple
      variables gets pretty tricky. The standard euclidean algorithm doesn't
      work for multivariate polynomials. Instead, we'd need to implement
      Grobner bases. I'm punting on this for the time being."""
      raise NotImplementedError

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
