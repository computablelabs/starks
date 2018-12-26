import unittest
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over

class TestNumberType(unittest.TestCase):
  """"
  Basic tests for classes in numbertype
  """

  def test_typecast(self):
    """Test typecasting operation"""
    mod3 = IntegersModP(3)
    Polynomial = polynomials_over(mod3)
    x = mod3(1)
    p = Polynomial([1,2])

    x+p
    p+x

