import unittest
import time
from starks.utils import mimc
from starks.merkle_tree import merkelize
from starks.merkle_tree import verify_branch
from starks.merkle_tree import mk_branch
from starks.merkle_tree import blake
from starks.merkle_tree import merkelize_polynomials
from starks.merkle_tree import unpack_merkle_leaf
from starks.air import Computation
from starks.air import get_computational_trace
from starks.fft import fft 
from starks.fri import prove_low_degree 
from starks.fri import verify_low_degree_proof 
from starks.utils import get_pseudorandom_indices
from starks.stark import get_power_cycle 
from starks.stark import mk_proof
from starks.stark import StarkParams 
from starks.stark import verify_proof
from starks.stark import construct_computation_polynomial 
from starks.stark import construct_constraint_polynomial 
from starks.stark import construct_remainder_polynomial 
from starks.stark import construct_boundary_polynomial 
#from starks.stark import construct_constants_polynomials
from starks.stark import get_pseudorandom_ks
from starks.stark import compute_pseudorandom_linear_combination_1d
from starks.stark import compute_pseudorandom_linear_combination
from starks.stark import compute_merkle_spot_checks
from starks.modp import IntegersModP
from starks.compression import bin_length
from starks.compression import compress_branches
from starks.compression import compress_fri
from starks.polynomial import polynomials_over
from starks.poly_utils import lagrange_interp_2
from starks.poly_utils import multivariates_over


class TestStark(unittest.TestCase):
  """
  Basic tests for Stark construction implementation. 
  """

  def test_stark_params(self):
    """Test generation of stark parameters."""
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    extension_factor = 8
    # Only tests that constructor works implicitly
    params = StarkParams(field, steps, modulus, extension_factor)

  def test_get_pseudorandom_indices(self):
    """
    Tests that pseudorandom indices are computed correctly.
    """
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    # state = [c_0, c_1, value]
    inp = [field(2), field(2), field(5)]
    width = 3
    steps = 512
    spot_check_security_factor = 80
    polysOver = multivariates_over(field, width).factory
    extension_factor = 8
    ## Factoring out computation
    # c_1*value**2 + c_0
    # Polys are: [X_1, X_2, X_1 + X_2*X_3**2]
    X_1 = polysOver({(1,0,0): field(1)})
    X_2 = polysOver({(0,1,0): field(1)})
    X_3 = polysOver({(0,0,1): field(1)})
    step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
    comp = Computation(field, width, inp, steps, step_polys, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    polys = [p_evaluations, d_evaluations, b_evaluations]
    mtree = merkelize_polynomials(width, polys)
    l_evaluations = compute_pseudorandom_linear_combination(
        comp, params, mtree, polys)
    l_mtree = merkelize(l_evaluations)

    indices = get_pseudorandom_indices(l_mtree[1],
        params.precision, count=spot_check_security_factor,
        exclude_multiples_of=params.extension_factor)
    assert len(indices) == spot_check_security_factor

  def test_compute_merkle_spot_checks(self):
    """
    Tests that merkle spot checks are constructed correctly.
    """
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    width = 3
    # state = [c_0, c_1, value]
    inp = [field(2), field(2), field(5)]
    steps = 512
    spot_check_security_factor = 80
    extension_factor = 8
    polysOver = multivariates_over(field, width).factory
    # c_1*value**2 + c_0
    # Polys are: [X_1, X_2, X_1 + X_2*X_3**2]
    X_1 = polysOver({(1,0,0): field(1)})
    X_2 = polysOver({(0,1,0): field(1)})
    X_3 = polysOver({(0,0,1): field(1)})
    step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
    comp = Computation(field, width, inp, steps, step_polys, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    polys = [p_evaluations, d_evaluations, b_evaluations]
    mtree = merkelize_polynomials(width, polys)
    l_evaluations = compute_pseudorandom_linear_combination(
        comp, params, mtree, polys)
    l_mtree = merkelize(l_evaluations)

    # Compute merkle spot checks
    branches = compute_merkle_spot_checks(mtree, l_mtree,
        comp, params, samples=spot_check_security_factor)

    # Each spot check returns 3 branches, m[pos], m[pos+extension_factor], l[pos] 
    assert len(branches) == 3*spot_check_security_factor

  def test_compute_pseudorandom_combination_1d(self):
    """
    Tests compute pseudorandom linear combination for 1-dimension
    """
    width = 3
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = IntegersModP(modulus)
    # state = [c_0, c_1, value]
    inp = [field(2), field(2), field(5)]
    polysOver = multivariates_over(field, width).factory
    ## Factoring out computation
    # c_1*value**2 + c_0
    # Polys are: [X_1, X_2, X_1 + X_2*X_3**2]
    X_1 = polysOver({(1,0,0): field(1)})
    X_2 = polysOver({(0,1,0): field(1)})
    X_3 = polysOver({(0,0,1): field(1)})
    step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
    comp = Computation(field, width, inp, steps, step_polys, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    polys = [p_evaluations, d_evaluations, b_evaluations]
    mtree = merkelize_polynomials(width, polys)
    l_evaluations_per_dim = compute_pseudorandom_linear_combination_1d(
        comp, params, mtree, polys)
    m_root = mtree[1]

    for i, pos in enumerate(range(params.precision)):
      p_of_x = p_evaluations[pos]
      next_pos = (pos + params.extension_factor) % params.precision
      p_of_g1x = p_evaluations[next_pos]
      d_of_x = d_evaluations[pos]
      b_of_x = b_evaluations[pos]

      x = params.G2**pos
      x_to_the_steps = x**comp.steps

      # Leaf node from l[pos]
      k1, k2, k3, k4 = get_pseudorandom_ks(mtree[1], 4)
      for dim in range(comp.width):
        l_evaluations_dim = l_evaluations_per_dim[dim]
        l_of_x_dim = l_evaluations_dim[pos]
        assert (l_of_x_dim - d_of_x[dim] - k1 * p_of_x[dim] - k2 * p_of_x[dim] * x_to_the_steps - k3 * b_of_x[dim] - k4 * b_of_x[dim] * x_to_the_steps) == 0

  def test_compute_pseudorandom_combination(self):
    """
    Tests compute pseudorandom linear combination 
    """
    width = 3
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = IntegersModP(modulus)
    inp = [field(2), field(2), field(5)]
    ## Factoring out computation
    # c_1*value**2 + c_0
    # Polys are: [X_1, X_2, X_1 + X_2*X_3**2]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0,0): field(1)})
    X_2 = polysOver({(0,1,0): field(1)})
    X_3 = polysOver({(0,0,1): field(1)})
    step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    polys = [p_evaluations, d_evaluations, b_evaluations]
    mtree = merkelize_polynomials(width, polys)
    l_evaluations_per_dim = compute_pseudorandom_linear_combination_1d(
        comp, params, mtree, polys)
    l_evaluations = compute_pseudorandom_linear_combination_1d(
        comp, params, mtree, polys)
    m_root = mtree[1]

    for i, pos in enumerate(range(params.precision)):
      p_of_x = p_evaluations[pos]
      next_pos = (pos + params.extension_factor) % params.precision
      p_of_g1x = p_evaluations[next_pos]
      d_of_x = d_evaluations[pos]
      b_of_x = b_evaluations[pos]

      x = params.G2**pos
      x_to_the_steps = x**comp.steps

      # Leaf node from l[pos]
      k1, k2, k3, k4 = get_pseudorandom_ks(mtree[1], 4)
      for dim in range(comp.width):
        l_evaluations_dim = l_evaluations_per_dim[dim]
        l_of_x_dim = l_evaluations_dim[pos]
        assert (l_of_x_dim - d_of_x[dim] - k1 * p_of_x[dim] - k2 * p_of_x[dim] * x_to_the_steps - k3 * b_of_x[dim] - k4 * b_of_x[dim] * x_to_the_steps) == 0


  def test_1d_end_to_end(self):
    """
    Tests stark end-to-end for 1d example 
    """
    width = 3
    steps = 512
    spot_check_security_factor = 80
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = IntegersModP(modulus)
    inp = [field(2), field(2), field(5)]
    ## Factoring out computation
    # c_1*value**2 + c_0
    # Polys are: [X_1, X_2, X_1 + X_2*X_3**2]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0,0): field(1)})
    X_2 = polysOver({(0,1,0): field(1)})
    X_3 = polysOver({(0,0,1): field(1)})
    step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    polys = [p_evaluations, d_evaluations, b_evaluations]
    mtree = merkelize_polynomials(width, polys)
    l_evaluations = compute_pseudorandom_linear_combination(
        comp, params, mtree, polys)
    l_mtree = merkelize(l_evaluations)
    m_root = mtree[1]
    l_root = l_mtree[1]
    fri_proof = prove_low_degree(
          l_evaluations,
          params.G2,
          # Manually filling in degree here
          steps * 3,
          modulus,
          exclude_multiples_of=extension_factor)

    assert verify_low_degree_proof(
        l_root,
        params.G2,
        fri_proof,
        # Manually filling in degree here
        steps * 3,
        modulus,
        exclude_multiples_of=extension_factor)

    # TODO(rbharath): Factor this into function?
    byte_list = [b'0x01', b'0x02', b'0x03', b'0x04']
    ks = [int.from_bytes(blake(mtree[1] + byte_list[ind]), 'big') for ind in range(4)]
    k1, k2, k3, k4 = ks

    branches = compute_merkle_spot_checks(mtree, l_mtree, comp, params)
    positions = get_pseudorandom_indices(l_mtree[1],
        params.precision, count=spot_check_security_factor,
        exclude_multiples_of=params.extension_factor)

    for i, pos in enumerate(positions):
      branch1 = mk_branch(mtree, pos)
      mbranch1 = verify_branch(m_root, pos, branch1)
      unpacked_leaf1 = unpack_merkle_leaf(mbranch1, comp.width, 3)
      next_pos = (pos + params.extension_factor) % params.precision
      branch2 = mk_branch(mtree, next_pos) 
      mbranch2 = verify_branch(m_root, next_pos, branch2)
      unpacked_leaf2 = unpack_merkle_leaf(mbranch2, comp.width, 3)
      # Leaf node from l[pos]
      l_of_x = verify_branch(l_root, pos, branches[i * 3 + 2],
          output_as_int=True)

      # Check that p_of_x was recovered correctly
      p_of_x = p_evaluations[pos]
      p_of_x_recovered = [field(p_of_x_dim) for p_of_x_dim in unpacked_leaf1[:comp.width]]
      assert len(p_of_x) == len(p_of_x_recovered)
      for dim in range(comp.width):
        assert p_of_x[dim] == p_of_x_recovered[dim]

      # Check that p_of_g1x was recovered correctly
      p_of_g1x = p_evaluations[next_pos]
      p_of_g1x_recovered = [field(p_of_g1x_dim) for p_of_g1x_dim in unpacked_leaf2[:comp.width]]
      assert len(p_of_g1x) == len(p_of_g1x_recovered)
      for dim in range(comp.width):
        assert p_of_g1x[dim] == p_of_g1x_recovered[dim]

      d_of_x = d_evaluations[pos]
      d_of_x_recovered = [field(d_of_x_dim) for d_of_x_dim in unpacked_leaf1[comp.width:2*comp.width]]
      assert len(d_of_x) == len(d_of_x_recovered)
      for dim in range(comp.width):
        assert d_of_x[dim] == d_of_x_recovered[dim]

      b_of_x = b_evaluations[pos]
      b_of_x_recovered = [field(b_of_x_dim) for b_of_x_dim in unpacked_leaf1[2*comp.width:]]
      assert len(b_of_x) == len(b_of_x_recovered)
      for dim in range(comp.width):
        assert b_of_x[dim] == b_of_x_recovered[dim]

      x = params.G2**pos
      x_to_the_steps = x**comp.steps
      zvalue = (x**comp.steps - 1)/(x - params.last_step_position)
      f_of_p_of_x = [comp.step_polys[i](p_of_x) for i in range(width)]
      f_of_p_of_x_recovered = [comp.step_polys[i](p_of_x) for i in range(width)]
      assert f_of_p_of_x == f_of_p_of_x_recovered
      assert (p_of_g1x[0] - f_of_p_of_x[0] - zvalue * d_of_x[0]) == 0
      assert (p_of_g1x_recovered[0] - f_of_p_of_x_recovered[0] - zvalue * d_of_x_recovered[0]) == 0

      # TODO(rbharath): How do I do this with multidimensional polynomial?
     # zeropoly2 = polysOver([-1, 1]) * polysOver([-params.last_step_position, 1])
     # for dim in range(comp.width):
     #   interpolant_dim = lagrange_interp_2(modulus, [1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
     #   assert (p_of_x[dim] - b_of_x[dim] * zeropoly2(x) - interpolant_dim(x)) == 0

  def test_higher_dim_trace(self):
    """
    Checks trace generation for multidimensional state.
    """
    width = 2
    steps = 5
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    trace, output = get_computational_trace(inp, steps,
        width, step_polys)
    assert list(trace[0]) == [0, 1]
    assert list(trace[1]) == [1, 1]
    assert list(trace[2]) == [1, 2]
    assert list(trace[3]) == [2, 3]
    assert list(trace[4]) == [3, 5]

  def test_higher_dim_computation_polynomial(self):
    """
    Tests construction of multidim computation polynomial
    """
    width = 2
    steps = 512
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    assert len(comp_poly_evals) == steps * extension_factor
    for cval in comp_poly_evals:
      assert isinstance(cval, list)
      assert len(cval) == width 


  def test_higher_dim_constraint_polynomial(self):
    """
    Tests construction of constraint polynomial.
    """
    width = 2
    steps = 512
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    constraint_evals = construct_constraint_polynomial(
        comp, params, comp_poly_evals)
    assert len(constraint_evals) == steps * extension_factor
    for cval in constraint_evals:
      assert isinstance(cval, list)
      assert len(cval) == width 

  def test_higher_dim_remainder_polynomial(self):
    """
    Basic tests of FRI generation for fibonacci stark
    """
    width = 2
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = IntegersModP(modulus)
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    inp = [field(0), field(1)]
    ## Factoring out computation
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    assert len(d_evaluations) == params.precision
    for ind, dval in enumerate(d_evaluations):
      assert isinstance(dval, list)
      assert len(dval) == width 
      for dim in range(width):
        assert isinstance(dval[dim], field)

  def test_higher_dim_boundary_polynomial(self):
    """
    Basic tests of FRI generation for fibonacci stark
    """
    width = 2
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    ## Factoring out computation
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)
    assert len(b_evaluations) == params.precision
    for ind, bval in enumerate(d_evaluations):
      assert isinstance(bval, list)
      assert len(bval) == width 
      for dim in range(width):
        assert isinstance(bval[dim], field)

  def test_higher_dim_fri(self):
    """
    Basic tests of FRI generation for fibonacci stark
    """
    width = 2
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    extension_factor = 8
    ## Factoring out computation
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    mtrees = merkelize_polynomials(width, [p_evaluations, d_evaluations, b_evaluations])

  def test_higher_dim_proof_verification(self):
    """
    Tests proof generation and verification for multidimensional state.
    """
    width = 2
    steps = 8
    constraint_degree = 4
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    # This is a place filler
    #constants = [[]] * steps
    extension_factor = 8
    #constraint_degree = 8
    #def step_fn(prev, constants):
    #  f_n_minus_1 = prev[0]
    #  f_n = prev[1]
    #  f_n_plus_1 = f_n + f_n_minus_1
    #  return [f_n, f_n_plus_1]
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_2, X_1 + X_2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    proof = mk_proof(comp, params)
    assert verify_proof(comp, params, proof)

  def test_computation_polynomial(self):
    """
    Tests construction of computation polynomial
    """
    dims = 1
    steps = 512
    extension_factor = 8
    constants = [[(i**7) ^ 42] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(5)]
    extension_factor = 8
    constraint_degree = 4
    def step_fn(state, constants):
      # 2value**2 + constant
      return [2*state[0]**2 + constants[0]]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    assert len(comp_poly_evals) == steps * extension_factor

  def test_constraint_polynomial(self):
    """
    Tests construction of constraint polynomial.
    """
    dims = 1
    steps = 512
    extension_factor = 8
    constants = [[(i**7) ^ 42] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(5)]
    extension_factor = 8
    constraint_degree = 4
    def step_fn(state, constants):
      # 2value**2 + constant
      return [2*state[0]**2 + constants[0]]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    constraint_evals = construct_constraint_polynomial(
        comp, params, comp_poly_evals)
    assert len(constraint_evals) == steps * extension_factor

  def test_compressed_stark(self):
    """Basic compressed stark test"""
    dims = 1
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(3)]
    extension_factor = 8
    constraint_degree = 4
    # Full STARK test
    constants = [[(i**7) ^ 42] for i in range(64)]
    constants = constants * (steps // 64)
    # Factoring out computation
    def step_fn(state, constants):
      return [state[0]**3 + constants[0]]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    proof = mk_proof(comp, params)
    m_root, l_root, branches, fri_proof = proof
    L1 = bin_length(compress_branches(branches))
    L2 = bin_length(compress_fri(fri_proof))
    print("Approx proof length: %d (branches), %d (FRI proof), %d (total)" %
          (L1, L2, L1 + L2))
    assert verify_proof(comp, params, proof)
                        
  def test_quadratic_stark(self):
    """
    Basic tests of quadratic stark generation
    """
    dims = 1
    steps = 512
    constraint_degree = 4
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    constants = [[(i**7) ^ 42] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(5)]
    extension_factor = 8

    # Factoring out computation
    def step_fn(state, constants):
      # 2value**2 + constant
      return [2*state[0]**2 + constants[0]]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    proof = mk_proof(comp, params)
    assert isinstance(proof, list)
    assert len(proof) == 4
    result = verify_proof(comp, params, proof)
    assert result


  def test_mimc_stark_verification(self):
    """
    Basic tests of MiMC stark verification.
    """
    dims = 1
    steps = 512
    constraint_degree = 4
    constants = [[(i**7) ^ 42] for i in range(64)]
    constants = constants * (steps // 64) 
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(5)]
    extension_factor = 8

    def step_fn(state, constants):
      return [state[0]**3 + constants[0]]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    proof = mk_proof(comp, params)
    result = verify_proof(comp, params, proof)
    assert result

  def test_affine_stark(self):
    """
    Basic tests of affine stark generation
    """
    dims = 1
    steps = 512
    constraint_degree = 4
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    constants = [[0]] * 512
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(5)]
    extension_factor = 8

    # Factoring out computation
    def step_fn(state, constants):
      return [3*state[0] + 4 + constants[0]]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    proof = mk_proof(comp, params)
    assert isinstance(proof, list)
    assert len(proof) == 4
    result = verify_proof(comp, params, proof)
    assert result

  def test_varying_quadratic_fri(self):
    """
    Basic tests of FRI generation for quadratic stark with varying coefficients
    """
    dims = 1
    steps = 512
    constraint_degree = 4
    constants = [[i, i] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = IntegersModP(modulus)
    inp = [field(5)]
    ## Factoring out computation
    def step_fn(state, constants):
      # c_1*value**2 + c_0
      value = state[0]
      return [constants[1]*value**2 + constants[0]]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)


    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    polys = [p_evaluations, d_evaluations, b_evaluations]
    mtree = merkelize_polynomials(dims, polys) 
    l_evaluations = compute_pseudorandom_linear_combination(
        comp, params, mtree, polys)
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
    dims = 1
    steps = 512
    constraint_degree = 4
    constants = [[i, i] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(5)]
    extension_factor = 8

    ## Factoring out computation
    def step_fn(value, constants):
      # c_1*value**2 + c_0
      return [constants[1]*value[0]**2 + constants[0]]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    proof = mk_proof(comp, params)
    assert isinstance(proof, list)
    assert len(proof) == 4
    result = verify_proof(comp, params, proof)
    assert result

  def test_varying_quintic_stark(self):
    """
    Basic tests of quintic stark generation
    """
    steps = 512
    dims = 1
    constraint_degree = 8
    constants = [[i]*6 for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(5)]
    extension_factor = 8

    ## Factoring out computation
    def step_fn(state, constants):
      # c_5*value**5 + c_4*value**4 + c_3*value**3 + c_2*value**2 + c_1*value**1 + c_0
      value = state[0]
      return [constants[5]*value**5 + constants[4]*value**4 + constants[3]*value**3 + constants[2]*value**2 + constants[1]*value + constants[0]]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, steps, modulus, extension_factor)
    proof = mk_proof(comp, params)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    trace, output = get_computational_trace(
        inp, steps, constants, step_fn)
    result = verify_proof(comp, params, proof)
    assert result
