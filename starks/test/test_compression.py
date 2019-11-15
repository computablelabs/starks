import unittest
from starks.compression import compress_fri
from starks.compression import decompress_fri
from starks.compression import bin_length
from starks.compression import compress_branches
from starks.fri import SmoothSubgroupFRI
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.utils import generate_Xi_s
from starks.air import AIR
from starks.stark import STARK

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
    root_of_unity = field(3)**((modulus-1)//8)

    fri = SmoothSubgroupFRI(field)
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree, modulus)
    compressed = compress_fri(proof)
    length = bin_length(compressed)
    print("bin_length: %d" % length)
    # TODO(rbharath): This is a lame test that checks length
    # of compressed proof is > 0. Need better unit test.
    assert length > 0

  #def test_compressed_stark(self):
  #  """Basic compressed stark test"""
  #  width = 2
  #  steps = 3
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  field = IntegersModP(modulus)
  #  inp = [field(2), field(5)]
  #  extension_factor = 8

  #  ## Factoring out computation
  #  [X_1, X_2] = generate_Xi_s(field, width)
  #  step_polys = [X_1, X_1 + X_2**3]
  #  comp = AIR(field, width, inp, steps, step_polys,
  #      extension_factor)
  #  stark = STARK(field, steps, extension_factor, width, step_polys)

  #  witness = comp.generate_witness()
  #  boundary = comp.generate_boundary_constraints()
  #  proof = stark.mk_proof(witness, boundary)
  #  m_root, l_root, branches, fri_proof = proof
  #  L1 = bin_length(compress_branches(branches))
  #  L2 = bin_length(compress_fri(fri_proof))
  #  print("Approx proof length: %d (branches), %d (FRI proof), %d (total)" %
  #        (L1, L2, L1 + L2))
  #  assert stark.verify_proof(proof, witness, boundary)
