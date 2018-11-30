import unittest
from starks.fft import fft

class TestFFT(unittest.TestCase):
  """
  Basic tests for fft implementation. 
  """
  def test_basic(self):
    """Basic test of fft."""
    modulus = 31 
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = [1, 2, 3, 4]
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = pow(3, (modulus-1)//6, modulus)
    evaluations = fft(poly, modulus, root_of_unity)
    assert len(evaluations) == 6

    # TODO(rbharath): This fails if root_of_unity = pow(3,
    # (modulus-1)//10, modulus) is a 10-th root of unity. Why?
    # It might have to do with an edge case bug in the code. 
