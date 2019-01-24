"""Test the algebraic placement and routing (APR) transformation."""
import unittest
from starks.modp import IntegersModP
from starks.finitefield import FiniteField 
from starks.air import Computation
from starks.apr import APR
from starks.poly_utils import multivariates_over

class TestAPR(unittest.TestCase):
  """Basic tests for APR classes."""

  def test_apr_constructor(self):
    """Test that the APR class can be initialized."""
    width = 2
    # Set the field small in tests since primitive polynomial generation is slow.
    p = 7 
    m = 4
    steps = 512
    extension_factor = 8
    field = FiniteField(p, m)
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    apr = APR(comp)


  def test_get_witness(self):
    """Test that a witness can be extracted from the APR instance.
    
    # TODO(rbharath): Make this test non-trivial to check the witness has
    # required properties.
    """
    width = 2
    # Set the field small in tests since primitive polynomial generation is slow.
    p = 7 
    m = 4
    steps = 4
    extension_factor = 8
    field = FiniteField(p, m)
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    apr = APR(comp)
    # TODO(rbharath): Uncomment this and reactivate it
    #witness = apr.generate_witness()
