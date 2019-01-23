import unittest
from starks.compression import compress_fri 
from starks.compression import decompress_fri 
from starks.compression import bin_length
from starks.fri import FRI
from starks.modp import IntegersModP
from starks.stark import StarkParams 
from starks.polynomial import polynomials_over

class TestCompression(unittest.TestCase):
  """
  Basic tests for compression/decompression of FRI proofs.
  """

  def test_compress_fri(self):
    """
    Basic tests of compression
    """
    degree = 4
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    polysOver = polynomials_over(field).factory
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = polysOver([val for val in range(degree)])
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = field(3)**((modulus-1)//6)
    #############################################
    from starks.fft import NonBinaryFFT
    fft_solver = NonBinaryFFT(field, root_of_unity)
    evals = fft_solver.fft(poly)
    assert 0 == 1
    #############################################
    
    fri = FRI(field, root_of_unity)
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree, modulus)
    compressed = compress_fri(proof)
    length = bin_length(compressed)
    print(compressed)
    print("bin_length: %d" % length)
    # TODO(rbharath): This is a lame test that checks length
    # of compressed proof is > 0. Need better unit test.
    assert length > 0
