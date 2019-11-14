"""Test the algebraic placement and routing (APR) transformation."""
import unittest
from starks.modp import IntegersModP
from starks.finitefield import FiniteField
from starks.air import AIR
from starks.apr import APR
from starks.poly_utils import multivariates_over
from starks.polynomial import polynomials_over

class TestAPR(unittest.TestCase):
  """Basic tests for APR classes."""

  def test_apr_constructor(self):
    """Test that the APR class can be initialized."""
    width = 2
    # Set the field small in tests since primitive polynomial generation is slow.
    p = 2
    m = 2
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    field = FiniteField(p, m)
    #m = 17
    ##x^17 + x^3 + 1 is primitive
    #coefficients = [Zp(0)] * 18
    #coefficients[0] = Zp(1)
    #coefficients[3] = Zp(1)
    #coefficients[17] = Zp(1)
    #poly = polysOver(coefficients)
    #field = FiniteField(p, m, polynomialModulus=poly)
    steps = 7
    extension_factor = 8
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2]
    air = AIR(field, width, inp, steps, step_polys,
        extension_factor)
    apr = APR(air)


  def test_get_witness(self):
    """Test that a witness can be extracted from the APR instance.

    # TODO(rbharath): Make this test non-trivial to check the witness has
    # required properties.
    """
    width = 2
    # Set the field small in tests since primitive polynomial generation is slow.
    p = 2
    m = 17
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    #x^17 + x^3 + 1 is primitive
    coefficients = [Zp(0)] * 18
    coefficients[0] = Zp(1)
    coefficients[3] = Zp(1)
    coefficients[17] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    steps = 3
    extension_factor = 8
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2]
    air = AIR(field, width, inp, steps, step_polys,
        extension_factor)
    apr = APR(air)
    # TODO(rbharath): Uncomment this and reactivate it
    witness = apr.generate_witness()
