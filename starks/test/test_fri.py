import unittest
from starks.fri import FRI
from starks.merkle_tree import merkelize
from starks.compression import bin_length
from starks.compression import compress_fri
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over


class TestFRI(unittest.TestCase):
  """
  Basic tests for FRI implementation. 
  """

  def test_basic_prove(self):
    """Test proof on low degree implementation"""
    degree = 4
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    polysOver = polynomials_over(field).factory
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = polysOver([val for val in range(degree)])

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = field(7)**((modulus-1)//8)

    # This is a low degree polynomial so we hit the special
    # case of the handler.
    fri = FRI(field, root_of_unity)
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree, modulus)
    # The proof is a list of length one, whose first entry is just the evaluations converted to bytes
    assert len(proof[0]) == 8 

  def test_high_degree_prove(self):
    """Tests proof generation on high degree polynomials"""
    steps = 512 
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    polysOver = polynomials_over(field).factory
    # Some round constants borrowed from MiMC
    poly = polysOver([field((i**7) ^ 42) for i in range(steps)])
    # Root of unity such that x^steps=1
    root_of_unity = field(7)**((modulus-1)//steps)
    # We're trying to prove this is a (steps-1)-degree
    # polnomial
    # degree = (steps-1) + 1 = steps
    fri = FRI(field)
    degree = steps
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree, modulus)
    # The proof recurses by dividing maxdeg_plus_1 by 4
    # So 512, 128, 32, 8. (The base case passes over to
    # special handler for degree 16 or less so these are all
    # recursions).
    assert len(proof) == 4
    for i, rec_proof in enumerate(proof):
      if i < 3:
        # Each subproof is [merkle_root, branches] for all but
        # base case.
        assert len(rec_proof) == 2
        assert len(rec_proof[1]) == 40
      else:
        # Here we trigger the base case.
        assert len(rec_proof) == 8

  def test_verify_low_degree_proof(self):
    """Verify a low degree proof"""
    dims = 1
    modulus = 31
    steps = 512
    degree = steps
    modulus = 2**256 - 2**32 * 351 + 1
    mod = IntegersModP(modulus)
    poly = [[mod((i**7) ^ 42)] for i in range(steps)]
    # Root of unity such that x^steps=1
    G = mod(7)**((modulus - 1) // steps)
    evaluations = fft(poly, modulus, G, dims=dims)
    # Unwrap the fft wrapping
    evaluations = [val[0] for val in evaluations]
    e_mtree = merkelize(evaluations)
    # This is a low degree polynomial so we hit the special
    # case of the handler.
    proof = prove_low_degree(evaluations, G, degree, modulus)

    root = e_mtree[1]
    verification = verify_low_degree_proof(root, G, proof, steps, modulus)

  def test_fri(self):
    """Pure FRI tests"""
    modulus = 2**256 - 2**32 * 351 + 1
    mod = IntegersModP(modulus)
    poly = [[mod(val)] for val in list(range(4096))]
    root_of_unity = mod(7)**((modulus - 1) // 16384)
    evaluations = fft(poly, modulus, root_of_unity)
    evaluations = [val[0] for val in evaluations]
    proof = prove_low_degree(evaluations, root_of_unity, 4096, modulus)
    print("Approx proof length: %d" % bin_length(compress_fri(proof)))
    assert verify_low_degree_proof(
        merkelize(evaluations)[1], root_of_unity, proof, 4096, modulus)

    try:
      fakedata = [
          x if pow(3, i, 4096) > 400 else 39 for x, i in enumerate(evaluations)
      ]
      proof2 = prove_low_degree(fakedata, root_of_unity, 4096, modulus)
      assert verify_low_degree_proof(
          merkelize(fakedata)[1], root_of_unity, proof, 4096, modulus)
      raise Exception("Fake data passed FRI")
    except:
      pass
    try:
      assert verify_low_degree_proof(
          merkelize(evaluations)[1], root_of_unity, proof, 2048, modulus)
      raise Exception("Fake data passed FRI")
    except:
      pass
