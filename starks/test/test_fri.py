import unittest
from starks.fri import prove_low_degree
from starks.fft import fft
from starks.poly_utils import PrimeField

class TestFRI(unittest.TestCase):
  """
  Basic tests for FRI implementation. 
  """
  def test_basic_prove(self):
    """Test proof on low degree implementation"""
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
    
    # This is a low degree polynomial so we hit the special
    # case of the handler.
    proof = prove_low_degree(evaluations, root_of_unity, degree, modulus)
    # The proof is a list of length one, whose first entry is just the evaluations converted to bytes
    assert len(proof[0]) == len(evaluations)

  def test_high_degree_prove(self):
    """Tests proof generation on high degree polynomials"""
    steps = 512 
    # Some round constants borrowed from MiMC 
    poly = [(i**7) ^ 42 for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)
    # Root of unity such that x^steps=1
    G = f.exp(7, (modulus - 1) // steps)
    evaluations = fft(poly, modulus, G)
    degree = steps

    # We're trying to prove this is a 512-degree polnomial
    proof = prove_low_degree(evaluations, G, degree, modulus)
    print("len(proof[0])")
    print(len(proof[0]))
    assert 0 == 1
