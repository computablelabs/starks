import unittest
from starks.fft import fft
from starks.fft import mul_polys
from starks.poly_utils import PrimeField
from starks.modp import IntegersModP

class TestFFT(unittest.TestCase):
  """
  Basic tests for fft implementation. 
  """
  def test_basic(self):
    """Basic test of fft."""
    modulus = 31 
    mod31 = IntegersModP(31)
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = [[mod31(val)] for val in range(4)] 
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = mod31(3)**((modulus-1)//6)
    evaluations = fft(poly, modulus, root_of_unity)
    assert len(evaluations) == 6

  def test_fft_output_type(self):
    """The output of FFT should be in the field if input is in field."""
    modulus = 31 
    mod31 = IntegersModP(31)
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = [[mod31(val)] for val in range(4)] 
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = mod31(3)**((modulus-1)//6)
    evaluations = fft(poly, modulus, root_of_unity)
    assert len(evaluations) == 6
    for val in evaluations:
      assert isinstance(val[0], mod31)

  def test_fft_multidim(self):
    """Test FFT of multidimensional signal."""
    steps = 512 
    modulus = 2**256 - 2**32 * 351 + 1
    mod = IntegersModP(modulus)
    # [1 + 2x + 3x^2 + 4 x^3 mod 31,
    #  1 + 2x + 3x^2 + 4 x^3 mod 31]
    poly = [[mod(i), mod(i)] for i in range(4)]
    # Root of unity such that x^512=1
    G = mod(7)**((modulus - 1) // steps)
    evaluations = fft(poly, modulus, G, dims=2)
    assert len(evaluations) == steps
    assert len(evaluations[0]) == 2

  def test_constants(self):
    """Test FFT handling of constants."""
    steps = 256 
    modulus = 2**256 - 2**32 * 351 + 1
    mod = IntegersModP(modulus)
    extension_factor = 8
    # precision = 512
    precision = steps * extension_factor
    # Root of unity such that x^512=1
    G2 = mod(7)**((modulus - 1) // precision)
    # Root of unity such that x^64=1
    G1 = G2**extension_factor
    round_constants = [(i**7) ^ 42 for i in range(64)]
    round_constants = [[mod(val)] for val in round_constants]
    skips2 = steps // len(round_constants)
    # Root of unity such that x^32 = 1
    root_of_unity = G1**skips2
    constants_mini_poly = fft(
        round_constants, modulus, root_of_unity, inv=True)
    assert len(constants_mini_poly) == len(round_constants)

  def test_mul_polys(self):
    """Test multiplication of polynomials."""
    modulus = 2**256 - 2**32 * 351 + 1
    mod = IntegersModP(modulus)
    # Root of unity such that x^512=1
    root_of_unity = mod(7)**((modulus - 1) // 512)
    # a = b = 1 + 2x + 3x^2 + 4 x^3 mod 31
    a = [mod(val) for val in range(4)] 
    b = [mod(val) for val in range(4)] 
    prod = mul_polys(a, b, root_of_unity)

