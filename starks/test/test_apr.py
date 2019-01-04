"""Test the algebraic placement and routing (APR) transformation."""
import unittest
from starks.modp import IntegersModP
from starks.air import Computation
from starks.apr import APR

class TestAPR(unittest.TestCase):
  """Basic tests for APR classes."""

  def test_apr_constructor(self):
    """Test that the APR class can be initialized."""
    dims = 2
    modulus = 2**256 - 2**32 * 351 + 1
    steps = 512
    constants = [[]] * steps
    extension_factor = 8
    field = IntegersModP(modulus)
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


