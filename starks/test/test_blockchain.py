import unittest
from starks.blockchain import Blockchain

class TestBlockchain(unittest.TestCase):
  """
  Basic tests for simple blockchain implementation. 
  """

  def test_construction(self):
    """Test construction of blockchain."""
    b = Blockchain()


