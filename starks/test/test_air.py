"""Test the construction of the algebraic intermediate representaiton."""

import unittest
from starks.air import AIR 
from starks.air import get_computational_trace
from starks.finitefield import FiniteField
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.poly_utils import multivariates_over
from starks.utils import generate_Xi_s


class TestAIR(unittest.TestCase):
  """
  Basic tests for AIR construction. 
  """

  def test_higher_dim_trace(self):
    """
    Checks trace generation for multidimensional state.
    """
    width = 2
    steps = 5
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    trace, output = get_computational_trace(inp, steps,
        width, step_polys)
    assert list(trace[0]) == [0, 1]
    assert list(trace[1]) == [1, 1]
    assert list(trace[2]) == [1, 2]
    assert list(trace[3]) == [2, 3]
    assert list(trace[4]) == [3, 5]

  def test_trace_extraction(self):
    """
    Tests construction of trace for a AIR 
    """
    width = 2
    steps = 512
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    [X_1, X_2] = generate_Xi_s(field, width)
    step_polys = [X_2, X_1 + X_2] 
    air = AIR(field, width, inp, steps, step_polys,
               extension_factor)
    assert len(air.computational_trace) == 512
    for state in air.computational_trace:
      assert len(state) == width 
      for dim in range(width):
        assert isinstance(state[dim], field)

