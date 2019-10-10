import pytest
from starks.modp import IntegersModP
from starks.multivariate_polynomial import multivariates_over
#from starks.utils import generate_Xi_s
#from starks.polynomial import polynomials_over
#from starks.finitefield import FiniteField

def test_basic_construction():
    """"Test construction."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory
    # This should equal y
    y_poly = multi({(0, 1, 0): 1})
  
    # Test construction of a constant polynomial
    zero_poly = multi(0)
    zero_poly_2 = multi({})
    assert zero_poly == zero_poly_2
#
#def test_degree():
#  """Test computation of multivariate degree"""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  n = 3
#  # Let's make polynomials in (Z/7)[x, y, z]
#  multi = multivariates_over(mod7, n).factory
#  # This should equal y
#  y_poly = multi({(0, 1, 0): 1})
#  assert y_poly.degree() == 1
#
#  # This should equal x + y + z
#  sum_poly = multi({(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1})
#  assert sum_poly.degree() == 1
#
#  # This should equal x^2 + y^2 + z^2
#  sq_sum_poly = multi({(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1})
#  assert sq_sum_poly.degree() == 2
#
#  # This should equal xy + yz + zx
#  multi_poly = multi({(1, 1, 0): 1, (0, 1, 1): 1, (1, 0, 1): 1})
#  assert multi_poly.degree() == 2
#
#def test_eq():
#  """Test multivariate poly equality."""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  n = 3
#  # Let's make polynomials in (Z/7)[x, y, z]
#  multi = multivariates_over(mod7, n).factory
#
#  # This should equal 0
#  zero_poly = multi({})
#  assert zero_poly == zero_poly
#
#  # This should equal y
#  y_poly = multi({(0, 1, 0): 1})
#  assert y_poly == y_poly
#
#  assert zero_poly != y_poly
#
#def test_neg():
#  """Tests negation of multivariate polynomials."""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  n = 3
#  # Let's make polynomials in (Z/7)[x, y, z]
#  multi = multivariates_over(mod7, n).factory
#
#  # This should equal 0
#  zero_poly = multi({})
#  assert zero_poly == -zero_poly
#
#  # This should equal y
#  y_poly = multi({(0, 1, 0): 1})
#  # This should equal -y
#  neg_y_poly = multi({(0, 1, 0): -1})
#  assert neg_y_poly == - y_poly
#
#def test_add():
#  """Test multivariate polynomial addition."""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  n = 3
#  # Let's make polynomials in (Z/7)[x, y, z]
#  multi = multivariates_over(mod7, n).factory
#
#  # This should equal 0
#  zero_poly = multi({})
#  assert zero_poly == zero_poly + zero_poly
#
#  # This should equal y
#  y_poly = multi({(0, 1, 0): 1})
#  assert y_poly == y_poly + zero_poly
#
#def test_mul():
#  """Test multiplication of multivariate polynomials."""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  n = 3
#  # Let's make polynomials in (Z/7)[x, y, z]
#  multi = multivariates_over(mod7, n).factory
#
#  # This should equal 0
#  zero_poly = multi({})
#  assert zero_poly == zero_poly * zero_poly
#
#  # This should equal y
#  y_poly = multi({(0, 1, 0): mod7(1)})
#  # This is y^2
#  y_sq_poly = multi({(0, 2, 0): mod7(1)})
#  assert y_sq_poly == y_poly * y_poly
#
#  # This should equal x + y
#  x_y_poly = multi({(1, 0, 0): mod7(1), (0, 1, 0): mod7(1)})
#  # This should equal x^2 + 2xy + y^2
#  sq_x_y_poly = multi({(2, 0, 0): mod7(1), (1, 1, 0): mod7(2), (0, 2, 0): mod7(1)})
#  prod = x_y_poly*x_y_poly
#  assert sq_x_y_poly == prod 
#
#def test_cross_mul():
#  """Test multiplication with different types."""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  n = 3
#  # Let's make polynomials in (Z/7)[x, y, z]
#  multi = multivariates_over(mod7, n).factory
#
#  three = mod7(3)
#  # This should equal y
#  y_poly = multi({(0, 1, 0): mod7(1)})
#  three_y_poly = multi({(0, 1, 0): mod7(3)})
#
#  assert three_y_poly == three * y_poly
#
#def test_exponentiation():
#  """Tests exponentiation of multivariate polynomials."""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  n = 3
#  # Let's make polynomials in (Z/7)[x, y, z]
#  multi = multivariates_over(mod7, n).factory
#
#  # This should equal 0
#  zero_poly = multi({})
#  assert zero_poly == zero_poly**2
#
#  # This should equal y
#  y_poly = multi({(0, 1, 0): mod7(1)})
#  # This is y^2
#  y_sq_poly = multi({(0, 2, 0): mod7(1)})
#  assert y_sq_poly == y_poly**2
#
#  # This should equal x + y
#  x_y_poly = multi({(1, 0, 0): mod7(1), (0, 1, 0): mod7(1)})
#  # This should equal x^2 + 2xy + y^2
#  sq_x_y_poly = multi({(2, 0, 0): mod7(1), (1, 1, 0): mod7(2), (0, 2, 0): mod7(1)})
#  assert sq_x_y_poly == x_y_poly**2
#
#def test_call():
#  """Test calling the multivariate polynomial like a function."""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  n = 3
#  # Let's make polynomials in (Z/7)[x, y, z]
#  multi = multivariates_over(mod7, n).factory
#  # This should equal 0
#  zero_poly = multi({})
#  # Evaluate with x=1, y=1, z=1. This should equal 0
#  assert zero_poly((1, 1, 1)) == 0
#
#  # This should equal y
#  y_poly = multi({(0, 1, 0): mod7(1)})
#  # Evaluate with x=1, y=1, z=1. This should equal 1 
#  assert y_poly((1, 1, 1)) == 1
#
#  # This should equal x^2 + 2xy + y^2
#  sq_x_y_poly = multi({(2, 0, 0): mod7(1), (1, 1, 0): mod7(2), (0, 2, 0): mod7(1)})
#  # Evaluate with x=1, y=1, z=1. This should equal 4 
#  assert sq_x_y_poly((1, 1, 1)) == 4
#
#def test_composition():
#  """Test that multivariate polynomials can compose properly."""
#  modulus = 2**256 - 2**32 * 351 + 1
#  field = IntegersModP(modulus)
#  width = 2
#  multi = multivariates_over(field, width).factory
#  [X_1, X_2] = generate_Xi_s(field, width)
#  step_polys = [X_2, X_1 + 2*X_2**2] 
#
#  # This should equal y + 1
#  state_polys = [X_1, X_2]
#  for step_poly in step_polys:
#    next_step = step_poly(state_polys)
#    # Since we start from the original polynomials the "next step" poly is
#    # just the transition poly.
#    assert next_step == step_poly
#
#def test_finitefield():
#  """Test multivariates over finite fields."""
#  # This finite field is of size 2^17
#  p = 2
#  m = 17
#  Zp = IntegersModP(p)
#  polysOver = polynomials_over(Zp)
#  #field = FiniteField(p, m)
#  #x^17 + x^3 + 1 is primitive 
#  coefficients = [Zp(0)] * 18
#  coefficients[0] = Zp(1)
#  coefficients[3] = Zp(1)
#  coefficients[17] = Zp(1)
#  poly = polysOver(coefficients)
#  field = FiniteField(p, m, polynomialModulus=poly)
#  width = 2
#  inp = [field(0), field(1)]
#  polysOver = multivariates_over(field, width).factory
#  [X_1, X_2] = generate_Xi_s(field, width)
#  step_poly = X_2
#  poly = step_poly([X_1, X_2])
#  assert poly == X_2
