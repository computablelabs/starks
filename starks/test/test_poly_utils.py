import unittest
from starks.poly_utils import PrimeField 

class TestPrimeField(unittest.TestCase):
  """
  Basic tests for PrimeField. 
  """

  def test_basic(self):
    """Basic test"""
    field7 = PrimeField(7)
    # 12 % 7 == 5
    out = field7.add(6, 6)
    assert out == 5

    # 6^-1 = 6
    out = field7.inv(6)
    assert out == 6

  def test_multi_inv(self):
    """Test of faster multiple inverse method."""
    field7 = PrimeField(7)
    # 6^-1 = 6
    outs = field7.multi_inv([6, 6, 6])
    assert outs == [6, 6, 6]

  def test_add_polys(self):
    """Test that addition of polynomials works."""
    field7 = PrimeField(7)
    added = field7.add_polys([1, 2], [1, 2])
    assert added == [2, 4]

