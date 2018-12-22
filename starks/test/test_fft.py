import unittest
from starks.fft import fft
from starks.poly_utils import PrimeField

class TestFFT(unittest.TestCase):
  """
  Basic tests for fft implementation. 
  """
  def test_basic(self):
    """Basic test of fft."""
    modulus = 31 
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = [[val] for val in range(4)] 
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = pow(3, (modulus-1)//6, modulus)
    evaluations = fft(poly, modulus, root_of_unity)
    assert len(evaluations) == 6

  def test_fft_multidim(self):
    """Test FFT of multidimensional signal."""
    steps = 512 
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)
    # [1 + 2x + 3x^2 + 4 x^3 mod 31,
    #  1 + 2x + 3x^2 + 4 x^3 mod 31]
    poly = [[1, 1], [2, 2], [3, 3], [4, 4]]
    # Root of unity such that x^512=1
    G = f.exp(7, (modulus - 1) // steps)
    evaluations = fft(poly, modulus, G, dims=2)
    assert len(evaluations) == steps
    assert len(evaluations[0]) == 2


  def test_constants(self):
    """Test FFT handling of constants."""
    steps = 256 
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)
    extension_factor = 8
    # precision = 512
    precision = steps * extension_factor
    # Root of unity such that x^512=1
    G2 = f.exp(7, (modulus - 1) // precision)
    # Root of unity such that x^64=1
    G1 = f.exp(G2, extension_factor)
    round_constants = [(i**7) ^ 42 for i in range(64)]
    round_constants = [[val] for val in round_constants]
    # skips2 = 2
    skips2 = steps // len(round_constants)
    # Root of unity such that x^32 = 1
    root_of_unity = f.exp(G1, skips2)
    constants_mini_poly = fft(
        round_constants, modulus, root_of_unity, inv=True)
    assert len(constants_mini_poly) == len(round_constants)
