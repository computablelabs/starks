import unittest
from starks.utils import mimc
from starks.utils import get_power_cycle 
from starks.modp import IntegersModP
from starks.utils import plus_one

class TestUtils(unittest.TestCase):
  """
  Basic tests for utils functions 
  """
  def test_mimc(self):
    """
    Basic tests of MiMC.
    """
    inp = 5
    steps = 3
    round_constants = [2, 7]
    val = mimc(inp, steps, round_constants)

  def test_get_power_cycle(self):
    """Basic test for power cycle."""
    modulus = 31 
    mod = IntegersModP(modulus)
    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    #root_of_unity = pow(3, (modulus-1)//6, modulus)
    root_of_unity = mod(3)**((modulus-1)//6)
    cycle = get_power_cycle(root_of_unity, modulus)
    assert cycle == [1, 26, 25, 30, 5, 6]

  def test_plus_one(self):
    """Test of simple typed function."""
    assert plus_one(4) == 5
