import unittest
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.finitefield import FiniteField
from starks.floatingpoint import FloatingPoint

class TestFloatingPoint(unittest.TestCase):
  """"
  Basic tests for floating point operations.
  """

  def test_addition(self):
    """Basic test of floating point addition."""

    p = 2
    m = 4
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    coefficients = [Zp(0)] * 5
    coefficients[0] = Zp(1)
    coefficients[1] = Zp(1)
    coefficients[4] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    floating_point = FloatingPoint(field)

    assert floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) == floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) + floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0)


  def test_subtraction(self):
    """Basic test of floating point subraction."""
    
    p = 2
    m = 4
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    coefficients = [Zp(0)] * 5
    coefficients[0] = Zp(1)
    coefficients[1] = Zp(1)
    coefficients[4] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    floating_point = FloatingPoint(field)

    assert floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) == floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) - floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0)

  def test_multiplication(self):
    """Basic test of floating point multiplication."""   

    p = 2
    m = 4
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    coefficients = [Zp(0)] * 5
    coefficients[0] = Zp(1)
    coefficients[1] = Zp(1)
    coefficients[4] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    floating_point = FloatingPoint(field)

    assert floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) == floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) * floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0)

  def test_division(self):
    """Basic test of floating point division."""
    
    p = 2
    m = 4
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    coefficients = [Zp(0)] * 5
    coefficients[0] = Zp(1)
    coefficients[1] = Zp(1)
    coefficients[4] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    floating_point = FloatingPoint(field)

    assert floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) == floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) / floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0)

  def test_exponentiation(self):
    """Basic test of floating point exponentiation."""
    

    p = 2
    m = 4
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    coefficients = [Zp(0)] * 5
    coefficients[0] = Zp(1)
    coefficients[1] = Zp(1)
    coefficients[4] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    floating_point = FloatingPoint(field)

    assert floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) == floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0) ** floating_point(field(polysOver([1,1,1])), field(polysOver([1,1,1])), 0, 0)


