import unittest
from starks.fri import prove_low_degree
from starks.fft import fft

class TestFRI(unittest.TestCase):
  """
  Basic tests for FRI implementation. 
  """
  def test_basic_prove(self):
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
    print(proof)
    print(len(proof))
    print(len(proof[0]))
    assert 0 == 1
