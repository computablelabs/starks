import unittest
from starks.compression import compress_fri 
from starks.compression import decompress_fri 
from starks.compression import bin_length
from starks.fri import prove_low_degree
from starks.fft import fft

class TestMiMC(unittest.TestCase):
  """
  Basic tests for compression/decompression of FRI proofs.
  """

  def test_compress_fri(self):
    """
    Basic tests of compression
    """
    degree = 4
    modulus = 31 
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = list(range(degree))
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = pow(3, (modulus-1)//6, modulus)
    evaluations = fft(poly, modulus, root_of_unity)
    assert len(evaluations) == 6
    
    proof = prove_low_degree(evaluations, root_of_unity, degree, modulus)
    compressed = compress_fri(proof)
    length = bin_length(compressed)
    print(compressed)
    print("bin_length: %d" % length)
    # TODO(rbharath): This is a lame test that checks length
    # of compressed proof is > 0. Need better unit test.
    assert length > 0
