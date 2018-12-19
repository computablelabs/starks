import unittest
import time
import numpy as np
from starks.utils import mimc
from starks.merkle_tree import merkelize
from starks.fft import fft 
from starks.fri import prove_low_degree 
from starks.fri import verify_low_degree_proof 
from starks.stark import get_power_cycle 
from starks.stark import mk_proof
from starks.stark import Computation
from starks.stark import StarkParams 
from starks.stark import verify_proof
from starks.stark import get_computational_trace
from starks.stark import construct_computation_polynomial 
from starks.stark import construct_constraint_polynomial 
from starks.stark import construct_remainder_polynomial 
from starks.stark import construct_boundary_polynomial 
from starks.stark import compute_pseudorandom_linear_combination
from starks.poly_utils import PrimeField
from starks.compression import bin_length
from starks.compression import compress_branches
from starks.compression import compress_fri


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
      return np.array([f_n, f_n_plus_1])
    trace, output = get_computational_trace(inp, steps,
        constants, fibonacci_step)
    assert list(trace[0]) == [0, 1]
    assert list(trace[1]) == [1, 1]
    assert list(trace[2]) == [1, 2]
    assert list(trace[3]) == [2, 3]
    assert list(trace[4]) == [3, 5]

  # TODO(rbharath): Fix this
  def test_higher_dim_computation_polynomial(self):
    """
    Tests construction of multidim computation polynomial
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    # This is a place filler
    constants = [[1] * steps]
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return np.array([f_n, f_n_plus_1])
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    comp = Computation(inp, steps, constants, step_fn)
    params = StarkParams(comp, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params, dims=dims)
    assert len(comp_poly_evals) == steps * extension_factor
    assert len(comp_poly_evals[0]) == dims

  def test_higher_dim_constraint_polynomial(self):
    """
    Tests construction of constraint polynomial.

    TODO(rbharath): This is failing
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    # This is a place filler
    constants = [[1] * steps]
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return np.array([f_n, f_n_plus_1])
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    comp = Computation(inp, steps, constants, step_fn)
    params = StarkParams(comp, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params, dims=dims)
    constraint_evals = construct_constraint_polynomial(
        comp, params, comp_poly_evals)
    assert len(constraint_evals) == steps * extension_factor
    for cval in constraint_evals:
      print("cval")
      print(cval)
      assert isinstance(cval, list)
      assert len(cval) == dims

  def test_higher_dim_remainder_polynomial(self):
    """
    Basic tests of FRI generation for fibonacci stark
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    # This is a place filler
    constants = [[1] * steps]
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    f = PrimeField(modulus)
    # This is a place filler
    constants = [[1] * steps]
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return np.array([f_n, f_n_plus_1])
    ## Factoring out computation
    comp = Computation(inp, steps, constants, step_fn)
    params = StarkParams(comp, modulus, extension_factor)


    p_evaluations = construct_computation_polynomial(
        comp, params, dims=dims)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations, dims=dims)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations, dims=dims)
    assert len(d_evaluations) == params.precision
    for ind, dval in enumerate(d_evaluations):
      assert len(dval) == dims
      for dim in range(dims):
        assert isinstance(dval[dim], int)

  def test_higher_dim_boundary_polynomial(self):
    """
    Basic tests of FRI generation for fibonacci stark
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    # This is a place filler
    constants = [[1] * steps]
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    f = PrimeField(modulus)
    # This is a place filler
    constants = [[1] * steps]
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return np.array([f_n, f_n_plus_1])
    ## Factoring out computation
    comp = Computation(inp, steps, constants, step_fn)
    params = StarkParams(comp, modulus, extension_factor)


    p_evaluations = construct_computation_polynomial(
        comp, params, dims=dims)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations, dims=dims)

  def test_higher_dim_fri(self):
    """
    Basic tests of FRI generation for fibonacci stark
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    # This is a place filler
    constants = [[1] * steps]
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    f = PrimeField(modulus)
    # This is a place filler
    constants = [[1] * steps]
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return np.array([f_n, f_n_plus_1])
    ## Factoring out computation
    comp = Computation(inp, steps, constants, step_fn)
    params = StarkParams(comp, modulus, extension_factor)


    p_evaluations = construct_computation_polynomial(
        comp, params, dims=dims)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations, dims=dims)

    mtrees = []
    for dim in range(dims):
      for pval, dval, bval in zip(p_evaluations, d_evaluations, b_evaluations):
        print("type(pval[dim]), type(dval[dim]), type(bval[dim])")
        print(type(pval[dim]), type(dval[dim]), type(bval[dim]))
        byte_val = pval[dim].to_bytes(32, 'big') + dval[dim].to_bytes(32, 'big') + bval[dim].to_bytes(
              32, 'big')
      dim_mtree = merkelize([
          pval[dim].to_bytes(32, 'big') + dval[dim].to_bytes(32, 'big') + bval[dim].to_bytes(
              32, 'big')
          for pval, dval, bval in zip(p_evaluations,
            d_evaluations, b_evaluations)
      ])
      mtrees.append(dim_mtree)

  ## TODO(rbharath): Fix this
  #def test_higher_dim_proof(self):
  #  """
  #  Tests proof generation for multidimensional state.

  #  TODO(rbharath): This test fails!!
  #  """
  #  dims = 2
  #  inp = [0, 1]
  #  steps = 8
  #  # This is a place filler
  #  constants = [[1] * steps]
  #  def fibonacci_step(f, prev, constants):
  #    f_n_minus_1 = prev[0]
  #    f_n = prev[1]
  #    f_n_plus_1 = f.add(f_n, f_n_minus_1)
  #    return np.array([f_n, f_n_plus_1])
  #  proof = mk_proof(inp, steps, constants, fibonacci_step,
  #      dims=dims)

  def test_computation_polynomial(self):
    """
    Tests construction of computation polynomial
    """
    inp = 5
    steps = 512
    extension_factor = 8
    constants = [[(i**7) ^ 42 for i in range(steps)]]
    modulus = 2**256 - 2**32 * 351 + 1
    def step_fn(f, value, constants):
      # 2value**2 + constant
      return f.add(f.mul(f.exp(value, 2), 2), constants[0])
    comp = Computation(inp, steps, constants, step_fn)
    params = StarkParams(comp, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    assert len(comp_poly_evals) == steps * extension_factor

  def test_constraint_polynomial(self):
    """
    Tests construction of constraint polynomial.
    """
    inp = 5
    steps = 512
    extension_factor = 8
    constants = [[(i**7) ^ 42 for i in range(steps)]]
    modulus = 2**256 - 2**32 * 351 + 1
    def step_fn(f, value, constants):
      # 2value**2 + constant
      return f.add(f.mul(f.exp(value, 2), 2), constants[0])
    comp = Computation(inp, steps, constants, step_fn)
    params = StarkParams(comp, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    constraint_evals = construct_constraint_polynomial(
        comp, params, comp_poly_evals)
    assert len(constraint_evals) == steps * extension_factor

  def test_compressed_stark(self):
    """Basic compressed stark test"""
    inp = 3
    steps = 512
    # Full STARK test
    round_constants = [(i**7) ^ 42 for i in range(64)]
    skips2 = steps // 64
    constants = round_constants * skips2
    constants = [constants]
    # Factoring out computation
    def mimc_step(f, inp, constants):
      return f.add(f.exp(inp, 3), constants[0])
    proof = mk_proof(inp, steps, constants, mimc_step)
    m_root, l_root, branches, fri_proof = proof
    L1 = bin_length(compress_branches(branches))
    L2 = bin_length(compress_fri(fri_proof))
    print("Approx proof length: %d (branches), %d (FRI proof), %d (total)" %
          (L1, L2, L1 + L2))
    trace, output = get_computational_trace(
        inp, steps, constants, mimc_step)
    assert verify_proof(inp, steps, constants, output, proof,
        mimc_step)
                        

  def test_mimc_stark(self):
    """
    Basic tests of MiMC Stark generation
    """
    inp = 5
    steps = 512
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
    extension_factor = 8
    f = PrimeField(modulus)
    ## Factoring out computation
    def step_fn(f, value, constants):
      # c_1*value**2 + c_0
      return f.add(f.mul(constants[1], f.exp(value, 2)), constants[0])
    comp = Computation(inp, steps, constants, step_fn)
    params = StarkParams(comp, modulus, extension_factor)


    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    mtree = merkelize([
        pval.to_bytes(32, 'big') + dval.to_bytes(32, 'big') + bval.to_bytes(
            32, 'big')
        for pval, dval, bval in zip(p_evaluations,
          d_evaluations, b_evaluations)
    ])
    l_evaluations = compute_pseudorandom_linear_combination(
        comp, params, mtree, d_evaluations, p_evaluations,
        b_evaluations)
    l_mtree = merkelize(l_evaluations)
    l_root = l_mtree[1]
    fri_proof = prove_low_degree(
          l_evaluations,
          params.G2,
          steps * constraint_degree,
          modulus,
          exclude_multiples_of=extension_factor)

    assert verify_low_degree_proof(
        l_root,
        params.G2,
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
