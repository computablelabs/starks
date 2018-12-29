"""Test the construction of the algebraic intermediate representaiton."""

import unittest
from starks.air import Computation
from starks.modp import IntegersModP


class TestAIR(unittest.TestCase):
  """
  Basic tests for AIR construction. 
  """

  def test_trace_extraction(self):
    """
    Tests construction of trace for a Computation 
    """
    dims = 2
    steps = 512
    constants = [[]] * steps
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    constraint_degree = 4
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f_n + f_n_minus_1
      return [f_n, f_n_plus_1]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    assert len(comp.computational_trace) == 512
    for state in comp.computational_trace:
        assert len(state) == dims
        for dim in range(dims):
          print("state[dim]")
          print(state[dim])
          assert isinstance(state[dim], field)
