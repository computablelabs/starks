import unittest
from starks.modp import IntegersModP
from starks.multivariate_polynomial import multivariates_over

class TestMultiVariatePolynomial(unittest.TestCase):
  """"
  Basic tests for finite field multivariate polynomial construction.
  """

  def test_basic_construction(self):
    """"Test construction."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory
    # This should equal y
    y__poly = multi({(0, 1, 0): 1})

