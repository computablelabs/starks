import unittest
from starks.utils import get_power_cycle 

class TestUtils(unittest.TestCase):
  """
  Basic tests for utils functions 
  """
  def test_get_power_cycle(self):
    """Basic test for power cycle."""
    modulus = 31 
    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = pow(3, (modulus-1)//6, modulus)
    cycle = get_power_cycle(root_of_unity, modulus)
    assert cycle == [1, 26, 25, 30, 5, 6]
