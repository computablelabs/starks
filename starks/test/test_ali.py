"""Test the Algebraic linking IOP (ALI) protocol."""

import unittest
from starks.air import Computation
from starks.apr import APR
from starks.ali import ALIProver
from starks.finitefield import FiniteField 

class TestALI(unittest.TestCase):
  """Basic tests for ALI classes."""

  def test_ali_prover_constructor(construction):
    """Tests that ALIProver can be initialized."""
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
    ali_prover = ALIProver(apr)
