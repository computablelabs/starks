import unittest
from starks.stark import mimc
from starks.stark import mk_proof
from starks.stark import verify_proof

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
    round_constants = [(i**7) ^ 42 for i in range(64)]
    modulus = 2**256 - 2**32 * 351 + 1
    # Factoring out computation
    def mimc_step(inp, i):
      return ((inp**3 + round_constants[i % len(round_constants)]) % modulus)

    def update_fn(f, value):
      """The step update function, in GF(2^n)"""
      return f.exp(value, 3)

    proof = mk_proof(inp, steps, round_constants, mimc_step, update_fn)
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
    scale_constants = [2 for _ in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    # Factoring out computation
    def quadratic_step(inp, i):
      return ((scale_constants[i]*inp**2 + round_constants[i]) % modulus)

    def update_fn(f, value, i):
      """The step update function, in GF(2^n)"""
      return f.mul(f.exp(value, 2), 2)

    def transition_constraint(value):
      return 2*value**2

    proof = mk_proof(inp, steps, round_constants, quadratic_step, update_fn)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    def computation(inp, steps, step_fn):
      for i in range(steps-1):
        inp = step_fn(inp, i)
      return inp
    output = computation(inp, steps, quadratic_step)
    result = verify_proof(inp, steps, round_constants,
                          output, proof, transition_constraint)
    assert result

  def test_mimc_stark_verification(self):
    """
    Basic tests of MiMC stark verification.
    """
    inp = 5
    LOGSTEPS = 9 
    steps = 2**LOGSTEPS
    round_constants = [(i**7) ^ 42 for i in range(64)]
    modulus = 2**256 - 2**32 * 351 + 1
    def mimc_step(inp, i):
      return ((inp**3 + round_constants[i % len(round_constants)]) % modulus)

    def update_fn(f, value):
      """The step update function, in GF(2^n)"""
      return f.exp(value, 3)
    proof = mk_proof(inp, steps, round_constants, mimc_step, update_fn)

    # The actual MiMC result
    output = mimc(inp, steps, round_constants)
    def transition_constraint(value):
      return value**3
    result = verify_proof(inp, steps, round_constants,
                          output, proof, transition_constraint)
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
    # Factoring out computation
    def affine_step(inp, i):
      return ((3*inp + 4) % modulus)

    def update_fn(f, value):
      """The step update function, in GF(2^n)"""
      return f.add(f.mul(3, value), 4)

    def transition_constraint(value):
      return ((3*value + 4) % modulus)

    proof = mk_proof(inp, steps, round_constants, affine_step, update_fn)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    def computation(inp, steps, step_fn):
      for i in range(steps-1):
        inp = step_fn(inp, i)
      return inp
    output = computation(inp, steps, affine_step)
    result = verify_proof(inp, steps, round_constants,
                          output, proof, transition_constraint)
    assert result
