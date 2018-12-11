import unittest
import time
from starks.stark import mk_proof
from starks.stark import verify_proof
from starks.poly_utils import PrimeField

# TODO(rbharath): Wait, does Vitalik's blog post claim that
# the verifier complexity is linear; Nah looks like t*log(t)
# is optimal. Verifier complexity is O(log**2(t)) which should
# be pretty small even for very large computations.
# NOTE(rbharath): These starks here are not zero-knowledge I
# think. Will need to be added onto library later.
def mimc(inp, steps, round_constants):
  """Compute a MIMC permutation for some number of steps"""
  modulus = 2**256 - 2**32 * 351 + 1
  start_time = time.time()
  for i in range(steps - 1):
    inp = (inp**3 + round_constants[i % len(round_constants)]) % modulus
  print("MIMC computed in %.4f sec" % (time.time() - start_time))
  return inp


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
    def mimc_step(f, inp, constant):
      return f.add(f.exp(inp, 3), constant)

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
    #scale_constants = [2 for _ in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)
    # Factoring out computation
    def quadratic_step(f, value, constant):
      # 2value**2 + constant
      return f.add(f.mul(f.exp(value, 2), 2), constant)

    proof = mk_proof(inp, steps, [round_constants], quadratic_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    def computation(inp, steps, constants, step_fn):
      round_constants = constants[0]
      for i in range(steps-1):
        inp = step_fn(f, inp, round_constants[i])
      return inp
    output = computation(inp, steps, [round_constants], quadratic_step)
    result = verify_proof(inp, steps, round_constants,
                          output, proof, quadratic_step)
    #assert 0 == 1
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
    def mimc_step(f, inp, constant):
      return f.add(f.exp(inp, 3), constant)

    proof = mk_proof(inp, steps, [round_constants], mimc_step)

    # The actual MiMC result
    output = mimc(inp, steps, round_constants)
    result = verify_proof(inp, steps, round_constants,
                          output, proof, mimc_step)
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
    def affine_step(f, value, constant):
      return f.add(f.add(f.mul(3, value), 4), constant)

    proof = mk_proof(inp, steps, [round_constants], affine_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    def computation(inp, steps, constants, step_fn):
      round_constants = constants[0]
      for i in range(steps-1):
        inp = step_fn(f, inp, round_constants[i])
      return inp
    output = computation(inp, steps, [round_constants], affine_step)
    result = verify_proof(inp, steps, round_constants,
                          output, proof, affine_step)
    assert result

  def test_varying_quadratic_stark(self):
    """
    Basic tests of quadratic stark with varying coefficients
    """
    inp = 5
    LOGSTEPS = 9 
    steps = 2**LOGSTEPS
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    round_constants = [(i**7) ^ 42 for i in range(steps)]
    scale_constants = [i for i in range(steps)]
    constants = [round_constants, scale_constants]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)
    # Factoring out computation
    def quadratic_step(f, value, constants):
      # c_1*value**2 + c_0 
      return f.add(f.mul(f.exp(value, constants[1]), 2), constants[0])

    proof = mk_proof(inp, steps, constants, quadratic_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    def computation(inp, steps, constants, step_fn):
      round_constants = constants[0]
      for i in range(steps-1):
        inp = step_fn(f, inp, [constants[0][i], constants[1][i]])
      return inp
    output = computation(inp, steps, constants, quadratic_step)
    result = verify_proof(inp, steps, round_constants,
                          output, proof, quadratic_step)
    #assert 0 == 1
    assert result
