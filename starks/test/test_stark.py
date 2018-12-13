import unittest
import time
from starks.utils import mimc
from starks.stark import mk_proof
from starks.stark import verify_proof
from starks.stark import get_computational_trace
from starks.poly_utils import PrimeField


class TestStark(unittest.TestCase):
  """
  Basic tests for Stark construction implementation. 
  """

  def test_mimc(self):
    """
    Basic tests of MiMC.
    """
    inp = 5
    steps = 3
    round_constants = [2, 7]
    val = mimc(inp, steps, round_constants)

  #def test_stark():
  #  """Basic stark test"""
  #  INPUT = 3
  #  LOGSTEPS = 13
  #  # Full STARK test
  #  constants = [(i**7) ^ 42 for i in range(64)]
  #  proof = mk_mimc_proof(INPUT, 2**LOGSTEPS, constants)
  #  m_root, l_root, branches, fri_proof = proof
  #  L1 = bin_length(compress_branches(branches))
  #  L2 = bin_length(compress_fri(fri_proof))
  #  print("Approx proof length: %d (branches), %d (FRI proof), %d (total)" %
  #        (L1, L2, L1 + L2))
  #  assert verify_mimc_proof(3, 2**LOGSTEPS, constants,
  #                           mimc(3, 2**LOGSTEPS, constants), proof)

  def test_mimc_stark(self):
    """
    Basic tests of MiMC Stark generation
    """
    inp = 5
    LOGSTEPS = 9
    steps = 2**LOGSTEPS
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    constants = [(i**7) ^ 42 for i in range(64)]
    skips2 = steps // 64
    round_constants = constants * skips2

    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)

    # Factoring out computation
    def mimc_step(f, inp, constants):
      return f.add(f.exp(inp, 3), constants[0])

    proof = mk_proof(inp, steps, [round_constants], mimc_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    # TODO(rbharath): Add more tests on these components

  def test_quadratic_stark(self):
    """
    Basic tests of quadratic stark generation
    """
    inp = 5
    LOGSTEPS = 9
    steps = 2**LOGSTEPS
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    round_constants = [(i**7) ^ 42 for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)

    # Factoring out computation
    def quadratic_step(f, value, constants):
      # 2value**2 + constant
      return f.add(f.mul(f.exp(value, 2), 2), constants[0])

    proof = mk_proof(inp, steps, [round_constants],
                     quadratic_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    trace, output = get_computational_trace(
        inp, steps, [round_constants], quadratic_step)
    result = verify_proof(inp, steps, [round_constants],
                          output, proof, quadratic_step)
    assert result

  def test_mimc_stark_verification(self):
    """
    Basic tests of MiMC stark verification.
    """
    inp = 5
    LOGSTEPS = 9
    steps = 2**LOGSTEPS
    constants = [(i**7) ^ 42 for i in range(64)]
    skips2 = steps // 64
    round_constants = constants * skips2

    modulus = 2**256 - 2**32 * 351 + 1

    def mimc_step(f, inp, constants):
      return f.add(f.exp(inp, 3), constants[0])

    proof = mk_proof(inp, steps, [round_constants], mimc_step)

    # The actual MiMC result
    output = mimc(inp, steps, round_constants)
    result = verify_proof(inp, steps, [round_constants], output, proof, mimc_step)
    assert result

  def test_affine_stark(self):
    """
    Basic tests of affine stark generation
    """
    inp = 5
    steps = 512
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    round_constants = [0 for i in range(512)]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)

    # Factoring out computation
    def affine_step(f, value, constants):
      return f.add(f.add(f.mul(3, value), 4), constants[0])

    proof = mk_proof(inp, steps, [round_constants], affine_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    _, output = get_computational_trace(inp, steps, [round_constants], affine_step)
    result = verify_proof(inp, steps, [round_constants], output, proof,
                          affine_step)
    assert result

  def test_varying_quadratic_stark(self):
    """
    Basic tests of quadratic stark with varying coefficients
    """
    inp = 5
    steps = 512
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    round_constants = [(i**7) ^ 42 for i in range(steps)]
    scale_constants = [3 for i in range(steps)]
    constants = [round_constants, scale_constants]
    #constants = [round_constants]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)

    ## Factoring out computation
    def quadratic_step(f, value, constants):
      # c_1*value**2 + c_0
      return f.add(f.mul(f.exp(value, constants[1]), 2), constants[0])
    #def quadratic_step(f, value, constants):
    #  # 2value**2 + constant
    #  return f.add(f.mul(f.exp(value, 2), 2), constants[0])

    proof = mk_proof(inp, steps, constants, quadratic_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    _, output = get_computational_trace(inp, steps, constants, quadratic_step)
    result = verify_proof(inp, steps, constants, output, proof, quadratic_step)
    assert result
