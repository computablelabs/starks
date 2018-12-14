import unittest
import time
from starks.utils import mimc
from starks.merkle_tree import merkelize
from starks.fft import fft 
from starks.fri import prove_low_degree 
from starks.fri import verify_low_degree_proof 
from starks.stark import get_power_cycle 
from starks.stark import mk_proof
from starks.stark import verify_proof
from starks.stark import get_computational_trace
from starks.stark import construct_constraint_polynomial 
from starks.stark import compute_remainder_polynomial 
from starks.stark import compute_boundary_polynomial 
from starks.stark import compute_pseudorandom_linear_combination
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

  def test_higher_dimensional_trace(self):
    """
    Checks trace generation for multidimensional state.
    """
    inp = [0, 1]
    steps = 5
    # This is a place filler
    constants = [[1] * steps]
    def fibonacci_step(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    trace, output = get_computational_trace(inp, steps,
        constants, fibonacci_step)
    assert trace[0] == [0, 1]
    assert trace[1] == [1, 1]
    assert trace[2] == [1, 2]
    assert trace[3] == [2, 3]
    assert trace[4] == [3, 5]

  def test_higher_dimensional_proof(self):
    """
    Tests proof generation for multidimensional state.
    """
    inp = [0, 1]
    steps = 8
    # This is a place filler
    constants = [[1] * steps]
    def fibonacci_step(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    proof = mk_proof(inp, steps, constants, fibonacci_step, dims=2)

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

  def test_cubic_stark(self):
    """
    Basic tests of cubic stark generation
    """
    inp = 5
    steps = 512
    round_constants = [i for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)

    # Factoring out computation
    def cubic_step(f, value, constants):
      # x**3 + 2value**2 + constant
      return f.add(f.exp(value, 3), f.add(f.mul(f.exp(value, 2), 2), constants[0]))

    proof = mk_proof(inp, steps, [round_constants],
                     cubic_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    trace, output = get_computational_trace(
        inp, steps, [round_constants], cubic_step)
    result = verify_proof(inp, steps, [round_constants],
                          output, proof, cubic_step)
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

  # TODO(rbharath): This doesn't work. Understand why...
  def test_varying_quadratic_fri(self):
    """
    Basic tests of FRI generation for quadratic stark with varying coefficients
    """
    inp = 5
    steps = 512
    constraint_degree = 2
    round_constants = [i for i in range(steps)]
    scale_constants = [i for i in range(steps)]
    constants = [round_constants, scale_constants]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)

    extension_factor = 8
    precision = steps * extension_factor
    G2 = f.exp(7, (modulus - 1) // precision)
    G1 = f.exp(G2, extension_factor)

    xs = get_power_cycle(G2, modulus)
    last_step_position = xs[(steps - 1) * extension_factor]

    ## Factoring out computation
    def step_fn(f, value, constants):
      # c_1*value**2 + c_0
      return f.add(f.mul(constants[1], f.exp(value, 2)), constants[0])

    computational_trace, output = get_computational_trace(inp,
        steps, constants, step_fn)
    computational_trace_polynomial = fft(
        computational_trace, modulus, G1, inv=True)
    p_evaluations = fft(computational_trace_polynomial,
        modulus, G2)
    c_of_p_evaluations = construct_constraint_polynomial(
        steps, constants, G1, G2, precision,
        p_evaluations, step_fn)
    d_evaluations = compute_remainder_polynomial(xs,
        precision, steps, last_step_position,
        c_of_p_evaluations)
    b_evaluations = compute_boundary_polynomial(xs,
        last_step_position, inp, output, p_evaluations)
    mtree = merkelize([
        pval.to_bytes(32, 'big') + dval.to_bytes(32, 'big') + bval.to_bytes(
            32, 'big')
        for pval, dval, bval in zip(p_evaluations,
          d_evaluations, b_evaluations)
    ])
    l_evaluations = compute_pseudorandom_linear_combination(mtree, G2, steps, precision, d_evaluations, p_evaluations, b_evaluations)
    l_mtree = merkelize(l_evaluations)
    l_root = l_mtree[1]
    fri_proof = prove_low_degree(
          l_evaluations,
          G2,
          steps * constraint_degree,
          modulus,
          exclude_multiples_of=extension_factor)


    assert verify_low_degree_proof(
        l_root,
        G2,
        fri_proof,
        steps * constraint_degree,
        modulus,
        exclude_multiples_of=extension_factor)

  def test_varying_quadratic_stark(self):
    """
    Basic tests of varying quadratic stark generation
    """
    inp = 5
    steps = 512
    round_constants = [i for i in range(steps)]
    scale_constants = [i for i in range(steps)]
    constants = [round_constants, scale_constants]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)

    ## Factoring out computation
    def quadratic_step(f, value, constants):
      # c_1*value**2 + c_0
      return f.add(f.mul(constants[1], f.exp(value, 2)), constants[0])

    proof = mk_proof(inp, steps, constants, quadratic_step)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    trace, output = get_computational_trace(
        inp, steps, constants, quadratic_step)
    result = verify_proof(inp, steps, constants,
                          output, proof, quadratic_step)
    assert result

  def test_varying_quintic_stark(self):
    """
    Basic tests of quintic stark generation
    """
    inp = 5
    steps = 512
    constraint_degree = 8
    zero_constants = [i for i in range(steps)]
    one_constants = [i for i in range(steps)]
    two_constants = [i for i in range(steps)]
    three_constants = [i for i in range(steps)]
    four_constants = [i for i in range(steps)]
    five_constants = [i for i in range(steps)]
    constants = [zero_constants, one_constants, two_constants, three_constants,
                 four_constants, five_constants]
    modulus = 2**256 - 2**32 * 351 + 1
    f = PrimeField(modulus)

    ## Factoring out computation
    def quintic_step(f, value, constants):
      # c_5*value**5 + c_4*value**4 + c_3*value**3 + c_2*value**2 + c_1*value**1 + c_0
      return f.add(f.mul(constants[5], f.exp(value, 5)),
          f.add(f.mul(constants[4], f.exp(value, 4)),
            f.add(f.mul(constants[3], f.exp(value, 3)),
              f.add(f.mul(constants[2], f.exp(value, 2)),
                f.add(f.mul(constants[1], f.exp(value, 1)), constants[0])))))

    proof = mk_proof(inp, steps, constants, quintic_step,
        constraint_degree=constraint_degree)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    trace, output = get_computational_trace(
        inp, steps, constants, quintic_step)
    result = verify_proof(inp, steps, constants,
                          output, proof, quintic_step,
                          constraint_degree=constraint_degree)
    assert result
