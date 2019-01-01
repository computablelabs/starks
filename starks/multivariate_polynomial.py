"""Classes and functions to aid representatin of multivariate polynomials.

Multivariate polynomials are used to represent transitions between adjacent
computation states.
"""
from __future__ import annotations
from typing import List 
from typing import Dict
from typing import Callable
from typing import Any
from starks.poly_utils import construct_multivariate_coefficients
from starks.numbertype import memoize
from starks.numbertype import Field
from starks.numbertype import FieldElement
from starks.numbertype import MultiVarPoly

# TODO(rbharath): How does the memoization code actually work?
@memoize
def multivariates_over(field: Field, num_vars: int) -> MultiVarPoly:
  """Create a multivariate polynomial.

  Let F be the field and n = num_vars. Then the polynomial ring we're
  constructing here is F[X_1,...,X_n].
  """

  class MultivariatePolynomial(MultiVarPoly):
    operatorPrecedence = 2

    # TODO(rbharath): Using Any here isn't optimal. cls is meant to be a field type. 
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
      elif isinstance(c, field):
        self.coefficients = {(0,)*num_vars: c}
      else:
        self.coefficients = c

    def __len__(self):
      return len(self.coefficients)

    def __repr__(self):
      if self.is_zero():
        return '0'

      def power_tuple_to_string(power_tup):
        return "".join(["X_%d^%d" % (i, power) in enumerate(power_tup)])

      return ' + '.join([
          '%s x^%d' % (coeff, power_tuple_to_string(power_tup)) if power_tup != (0,)*num_vars else '%s' % coeff 
          for power_tup, coeff in enumerate(self.coefficients)
      ])

    def degree(self):
      """TODO(rbharath): Computing the degree is a little tricky."""

      raise NotImplementedError

    def __sub__(self, other):
      return self + (-other)

    def __iter__(self):
      return iter(self.coefficients)

    def __neg__(self):
      return Polynomial({(power_tup, -coeff) for (power_tup, coeff) in self})


  MultivariatePolynomial.field = field
  MultivariatePolynomial.num_vars = num_vars
  MultivariatePolynomial.__name__ = "".join(["(%s)" % field.__name__, "[", ",".join(["X_%d" % i for i in range(num_vars)]), "]"])
  return MultivariatePolynomial
