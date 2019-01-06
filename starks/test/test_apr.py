"""Test the algebraic placement and routing (APR) transformation."""
import unittest
from starks.modp import IntegersModP
from starks.finitefield import FiniteField 
from starks.air import Computation
from starks.apr import APR

class TestAPR(unittest.TestCase):
  """Basic tests for APR classes."""

  def test_apr_constructor(self):
    """Test that the APR class can be initialized."""
    dims = 2
    # Set the field small in tests since primitive polynomial generation is slow.
    p = 7 
    m = 4
    steps = 512
    constants = [[]] * steps
    extension_factor = 8
    field = FiniteField(p, m)
    inp = [field(0), field(1)]
    constraint_degree = 4
    def step_fn(prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f_n + f_n_minus_1
      return [f_n, f_n_plus_1]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    apr = APR(comp)


  def test_get_witness(self):
    """Test that a witness can be extracted from the APR instance.
    
    # TODO(rbharath): Make this test non-trivial to check the witness has
    # required properties.
    """
    dims = 2
    # Set the field small in tests since primitive polynomial generation is slow.
    p = 7 
    m = 4
    steps = 4
    constants = [[]] * steps
    extension_factor = 8
    field = FiniteField(p, m)
    inp = [field(0), field(1)]
    constraint_degree = 4
    def step_fn(prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f_n + f_n_minus_1
      return [f_n, f_n_plus_1]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    apr = APR(comp)
    witness = apr.generate_witness()
