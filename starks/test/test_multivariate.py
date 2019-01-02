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
    y_poly = multi({(0, 1, 0): 1})

  def test_degree(self):
    """Test computation of multivariate degree"""
    modulus = 7
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory
    # This should equal y
    y_poly = multi({(0, 1, 0): 1})
    assert y_poly.degree() == 1

    # This should equal x + y + z
    sum_poly = multi({(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1})
    assert sum_poly.degree() == 1

    # This should equal x^2 + y^2 + z^2
    sq_sum_poly = multi({(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1})
    assert sq_sum_poly.degree() == 2

    # This should equal xy + yz + zx
    multi_poly = multi({(1, 1, 0): 1, (0, 1, 1): 1, (1, 0, 1): 1})
    assert multi_poly.degree() == 2

  def test_eq(self):
    """Test multivariate poly equality."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory

    # This should equal 0
    zero_poly = multi({})
    assert zero_poly == zero_poly

    # This should equal y
    y_poly = multi({(0, 1, 0): 1})
    assert y_poly == y_poly

    assert zero_poly != y_poly

  def test_add(self):
    """Test multivariate polynomial addition."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory

    # This should equal 0
    zero_poly = multi({})
    assert zero_poly == zero_poly + zero_poly

    # This should equal y
    y_poly = multi({(0, 1, 0): 1})
    assert y_poly == y_poly + zero_poly

  def test_mul(self):
    """Test multiplication of multivariate polynomials."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory

    # This should equal 0
    zero_poly = multi({})
    assert zero_poly == zero_poly * zero_poly

    # This should equal y
    y_poly = multi({(0, 1, 0): mod7(1)})
    # This is y^2
    y_sq_poly = multi({(0, 2, 0): mod7(1)})
    assert y_sq_poly == y_poly * y_poly

    # This should equal x + y
    x_y_poly = multi({(1, 0, 0): mod7(1), (0, 1, 0): mod7(1)})
    # This should equal x^2 + 2xy + y^2
    sq_x_y_poly = multi({(2, 0, 0): mod7(1), (1, 1, 0): mod7(2), (0, 2, 0): mod7(1)})
    prod = x_y_poly*x_y_poly
    assert sq_x_y_poly == prod 

  def test_exponentiation(self):
    """Tests exponentiation of multivariate polynomials."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory

    # This should equal 0
    zero_poly = multi({})
    assert zero_poly == zero_poly**2

    # This should equal y
    y_poly = multi({(0, 1, 0): mod7(1)})
    # This is y^2
    y_sq_poly = multi({(0, 2, 0): mod7(1)})
    assert y_sq_poly == y_poly**2

    # This should equal x + y
    x_y_poly = multi({(1, 0, 0): mod7(1), (0, 1, 0): mod7(1)})
    # This should equal x^2 + 2xy + y^2
    sq_x_y_poly = multi({(2, 0, 0): mod7(1), (1, 1, 0): mod7(2), (0, 2, 0): mod7(1)})
    assert sq_x_y_poly == x_y_poly**2
