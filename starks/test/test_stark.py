import unittest
import time
from starks.utils import mimc
from starks.merkle_tree import merkelize
from starks.merkle_tree import verify_branch
from starks.merkle_tree import mk_branch
from starks.merkle_tree import blake
from starks.air import Computation
from starks.air import get_computational_trace
from starks.fft import fft 
from starks.fri import prove_low_degree 
from starks.fri import verify_low_degree_proof 
from starks.utils import get_pseudorandom_indices
from starks.stark import get_power_cycle 
from starks.stark import unpack_merkle_leaf
from starks.stark import mk_proof
from starks.stark import StarkParams 
from starks.stark import verify_proof
from starks.stark import construct_computation_polynomial 
from starks.stark import construct_constraint_polynomial 
from starks.stark import construct_remainder_polynomial 
from starks.stark import construct_boundary_polynomial 
from starks.stark import construct_constants_polynomials
from starks.stark import get_pseudorandom_ks
from starks.stark import compute_pseudorandom_linear_combination_1d
from starks.stark import compute_pseudorandom_linear_combination
from starks.stark import merkelize_polynomials 
from starks.stark import compute_merkle_spot_checks
#from starks.poly_utils import PrimeField
from starks.modp import IntegersModP
from starks.compression import bin_length
from starks.compression import compress_branches
from starks.compression import compress_fri


class TestStark(unittest.TestCase):
  """
  Basic tests for Stark construction implementation. 
  """

  def test_unpack_merkle_leaf(self):
    """
    Tests that merkle leaf unpacking works correctly.
    """
    leaf_parts = []
    dims = 2
    num_polys = 3
    count = 0
    for poly_ind in range(num_polys):
      for dim in range(dims):
        count_bytes = count.to_bytes(32, 'big')
        leaf_parts.append(count_bytes)
        count += 1
    leaf = b''.join(leaf_parts)
    vals = unpack_merkle_leaf(leaf, dims, num_polys)
    assert len(vals) == len(leaf_parts)
    for ind, val in enumerate(vals):
      assert val == leaf_parts[ind]

  def test_get_pseudorandom_indices(self):
    """
    Tests that pseudorandom indices are computed correctly.
    """
    dims = 1
    inp = [5]
    steps = 512
    constraint_degree = 4
    spot_check_security_factor = 80
    constants = [[i, i] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    extension_factor = 8
    ## Factoring out computation
    def step_fn(f, state, constants):
      # c_1*value**2 + c_0
      #return [f.add(f.mul(constants[1], f.exp(value, 2)), constants[0])]
      return [constants[1]*state[0]**2 + constants[0]]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)

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

    # TODO(rbharath): I'm copy-pasting test setup here. There
    # ought to be a simpler way to do this.
    # TODO(rbharath): Can simplify this test with just a
    # random seed.
    indices = get_pseudorandom_indices(l_mtree[1],
        params.precision, count=spot_check_security_factor,
        exclude_multiples_of=params.extension_factor)
    assert len(indices) == spot_check_security_factor

  def test_compute_merkle_spot_checks(self):
    """
    Tests that merkle spot checks are constructed correctly.
    """
    dims = 1
    inp = [5]
    steps = 512
    constraint_degree = 4
    spot_check_security_factor = 80
    constants = [[i, i] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = PrimeField(modulus)
    ## Factoring out computation
    def step_fn(f, value, constants):
      # c_1*value**2 + c_0
      return [f.add(f.mul(constants[1], f.exp(value[0], 2)), constants[0])]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)

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

    # Compute merkle spot checks
    branches = compute_merkle_spot_checks(mtree, l_mtree,
        comp, params, samples=spot_check_security_factor)

    # Each spot check returns 3 branches, m[pos], m[pos+extension_factor], l[pos] 
    assert len(branches) == 3*spot_check_security_factor

  def test_compute_pseudorandom_combination_1d(self):
    """
    Tests compute pseudorandom linear combination for 1-dimension
    """
    inp = [5]
    dims = 1
    steps = 512
    constants = [[i, i] for i in range(512)]
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    constraint_degree = 4
    field = PrimeField(modulus)
    ## Factoring out computation
    def step_fn(f, state, constants):
      # c_1*value**2 + c_0
      return [f.add(f.mul(constants[1], f.exp(state[0], 2)), constants[0])]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)

    constants_extensions, constants_polynomials = \
        construct_constants_polynomials(comp, params)
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
    l_evaluations_per_dim = compute_pseudorandom_linear_combination_1d(
        comp, params, mtree, polys)
    m_root = mtree[1]

    for i, pos in enumerate(range(params.precision)):
      p_of_x = p_evaluations[pos]
      next_pos = (pos + params.extension_factor) % params.precision
      p_of_g1x = p_evaluations[next_pos]
      d_of_x = d_evaluations[pos]
      b_of_x = b_evaluations[pos]

      x = field.exp(params.G2, pos)
      x_to_the_steps = field.exp(x, comp.steps)

      # Leaf node from l[pos]
      k1, k2, k3, k4 = get_pseudorandom_ks(mtree[1], 4)
      for dim in range(comp.dims):
        l_evaluations_dim = l_evaluations_per_dim[dim]
        l_of_x_dim = l_evaluations_dim[pos]
        assert (l_of_x_dim - d_of_x[dim] - k1 * p_of_x[dim] - k2 * p_of_x[dim] * x_to_the_steps - k3 * b_of_x[dim] - k4 * b_of_x[dim] * x_to_the_steps) % modulus == 0

  def test_compute_pseudorandom_combination(self):
    """
    Tests compute pseudorandom linear combination 
    """
    inp = [5]
    dims = 1
    steps = 512
    constants = [[i, i] for i in range(512)]
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    constraint_degree = 4
    field = PrimeField(modulus)
    ## Factoring out computation
    def step_fn(f, value, constants):
      # c_1*value**2 + c_0
      return [f.add(f.mul(constants[1], f.exp(value[0], 2)), constants[0])]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)

    constants_extensions, constants_polynomials = \
        construct_constants_polynomials(comp, params)
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

      x = field.exp(params.G2, pos)
      x_to_the_steps = field.exp(x, comp.steps)

      # Leaf node from l[pos]
      k1, k2, k3, k4 = get_pseudorandom_ks(mtree[1], 4)
      for dim in range(comp.dims):
        l_evaluations_dim = l_evaluations_per_dim[dim]
        l_of_x_dim = l_evaluations_dim[pos]
        assert (l_of_x_dim - d_of_x[dim] - k1 * p_of_x[dim] - k2 * p_of_x[dim] * x_to_the_steps - k3 * b_of_x[dim] - k4 * b_of_x[dim] * x_to_the_steps) % modulus == 0


  def test_1d_end_to_end(self):
    """
    Tests stark end-to-end for 1d example 
    """
    dims = 1
    inp = [5]
    steps = 512
    constraint_degree = 4
    constants = [[i, i] for i in range(steps)]
    spot_check_security_factor = 80
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = PrimeField(modulus)
    ## Factoring out computation
    def step_fn(f, value, constants):
      # c_1*value**2 + c_0
      return [f.add(f.mul(constants[1], f.exp(value[0], 2)), constants[0])]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)

    constants_extensions, constants_polynomials = \
        construct_constants_polynomials(comp, params)
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
    m_root = mtree[1]
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

    # TODO(rbharath): Factor this into function?
    byte_list = [b'0x01', b'0x02', b'0x03', b'0x04']
    ks = [int.from_bytes(blake(mtree[1] + byte_list[ind]), 'big') for ind in range(4)]
    k1, k2, k3, k4 = ks

    branches = compute_merkle_spot_checks(mtree, l_mtree, comp, params)
    #branches = []
    positions = get_pseudorandom_indices(l_mtree[1],
        params.precision, count=spot_check_security_factor,
        exclude_multiples_of=params.extension_factor)

    for i, pos in enumerate(positions):
      branch1 = mk_branch(mtree, pos)
      mbranch1 = verify_branch(m_root, pos, branch1)
      unpacked_leaf1 = unpack_merkle_leaf(mbranch1, comp.dims, 3)
      next_pos = (pos + params.extension_factor) % params.precision
      branch2 = mk_branch(mtree, next_pos) 
      mbranch2 = verify_branch(m_root, next_pos, branch2)
      unpacked_leaf2 = unpack_merkle_leaf(mbranch2, comp.dims, 3)
      # Leaf node from l[pos]
      l_of_x = verify_branch(l_root, pos, branches[i * 3 + 2],
          output_as_int=True)

      # Check that p_of_x was recovered correctly
      p_of_x = p_evaluations[pos]
      p_of_x_recovered = [int.from_bytes(p_of_x_dim, 'big') for p_of_x_dim in unpacked_leaf1[:comp.dims]]
      assert len(p_of_x) == len(p_of_x_recovered)
      for dim in range(comp.dims):
        assert p_of_x[dim] == p_of_x_recovered[dim]

      # Check that p_of_g1x was recovered correctly
      # TODO(rbharath): This is breaking down!!
      p_of_g1x = p_evaluations[next_pos]
      p_of_g1x_recovered = [int.from_bytes(p_of_g1x_dim, 'big') for p_of_g1x_dim in unpacked_leaf2[:comp.dims]]
      assert len(p_of_g1x) == len(p_of_g1x_recovered)
      for dim in range(comp.dims):
        assert p_of_g1x[dim] == p_of_g1x_recovered[dim]

      d_of_x = d_evaluations[pos]
      d_of_x_recovered = [int.from_bytes(d_of_x_dim, 'big') for d_of_x_dim in unpacked_leaf1[comp.dims:2*comp.dims]]
      assert len(d_of_x) == len(d_of_x_recovered)
      for dim in range(comp.dims):
        assert d_of_x[dim] == d_of_x_recovered[dim]

      b_of_x = b_evaluations[pos]
      b_of_x_recovered = [int.from_bytes(b_of_x_dim, 'big') for b_of_x_dim in unpacked_leaf1[2*comp.dims:]]
      assert len(b_of_x) == len(b_of_x_recovered)
      for dim in range(comp.dims):
        assert b_of_x[dim] == b_of_x_recovered[dim]

      x = field.exp(params.G2, pos)
      x_to_the_steps = field.exp(x, comp.steps)
      zvalue = field.div(field.exp(x, comp.steps) - 1, x - params.last_step_position)
      k_of_xs = []
      for constants_mini_polynomial in constants_polynomials:
        # This is unwrapping the polynomial
        constants_mini_polynomial = [val[0] for val in constants_mini_polynomial]
        k_of_x = field.eval_poly_at(constants_mini_polynomial, x)
        k_of_xs.append(k_of_x)
      f_of_p_of_x = comp.step_fn(field, p_of_x, k_of_xs)
      f_of_p_of_x_recovered = comp.step_fn(field, p_of_x, k_of_xs)
      assert f_of_p_of_x == f_of_p_of_x_recovered
      assert (p_of_g1x[0] - f_of_p_of_x[0] - zvalue * d_of_x[0]) % modulus == 0
      assert (p_of_g1x_recovered[0] - f_of_p_of_x_recovered[0] - zvalue * d_of_x_recovered[0]) % modulus == 0

      zeropoly2 = field.mul_polys([-1, 1], [-params.last_step_position, 1])
      for dim in range(comp.dims):
        interpolant_dim = field.lagrange_interp_2([1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
        assert (p_of_x[dim] - b_of_x[dim] * field.eval_poly_at(zeropoly2, x) - field.eval_poly_at(interpolant_dim, x)) % modulus == 0

  def test_higher_dim_trace(self):
    """
    Checks trace generation for multidimensional state.
    """
    inp = [0, 1]
    steps = 5
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    # This is a place filler
    constants = [[]] * steps
    def fibonacci_step(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    trace, output = get_computational_trace(field, inp, steps,
        constants, fibonacci_step)
    assert list(trace[0]) == [0, 1]
    assert list(trace[1]) == [1, 1]
    assert list(trace[2]) == [1, 2]
    assert list(trace[3]) == [2, 3]
    assert list(trace[4]) == [3, 5]

  def test_higher_dim_computation_polynomial(self):
    """
    Tests construction of multidim computation polynomial
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    constants = [[]] * steps
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    constraint_degree = 4
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    assert len(comp_poly_evals) == steps * extension_factor
    for cval in comp_poly_evals:
      assert isinstance(cval, list)
      assert len(cval) == dims

  def test_higher_dim_constants_polynomial(self):
    """
    Tests construction of multidim constants polynomial
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    constants = [[1]] * steps
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    constraint_degree = 4
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    extensions, polys = construct_constants_polynomials(
        comp, params)
    assert len(extensions) == 1
    assert len(polys) == 1
    assert len(extensions[0]) == params.precision
    assert len(polys[0]) == steps 


  def test_higher_dim_constraint_polynomial(self):
    """
    Tests construction of constraint polynomial.
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    # This is a place filler
    constants = [[]] * steps
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    constraint_degree = 4
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    constraint_evals = construct_constraint_polynomial(
        comp, params, comp_poly_evals)
    assert len(constraint_evals) == steps * extension_factor
    for cval in constraint_evals:
      assert isinstance(cval, list)
      assert len(cval) == dims

  def test_higher_dim_remainder_polynomial(self):
    """
    Basic tests of FRI generation for fibonacci stark
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    constants = [[]] * steps
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = PrimeField(modulus)
    constraint_degree = 4
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    ## Factoring out computation
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    assert len(d_evaluations) == params.precision
    for ind, dval in enumerate(d_evaluations):
      assert isinstance(dval, list)
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
    constants = [[]] * steps
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = PrimeField(modulus)
    constraint_degree = 4
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    ## Factoring out computation
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)

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
      assert len(bval) == dims
      for dim in range(dims):
        assert isinstance(bval[dim], int)

  def test_higher_dim_fri(self):
    """
    Basic tests of FRI generation for fibonacci stark
    """
    dims = 2
    inp = [0, 1]
    steps = 512
    constants = [[]] * steps
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = PrimeField(modulus)
    extension_factor = 8
    constraint_degree = 4
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    ## Factoring out computation
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)

    p_evaluations = construct_computation_polynomial(
        comp, params)
    c_of_p_evaluations = construct_constraint_polynomial(
        comp, params, p_evaluations)
    d_evaluations = construct_remainder_polynomial(
        comp, params, c_of_p_evaluations)
    b_evaluations = construct_boundary_polynomial(
        comp, params, p_evaluations)

    mtrees = merkelize_polynomials(dims, [p_evaluations, d_evaluations, b_evaluations])

  def test_higher_dim_proof_verification(self):
    """
    Tests proof generation and verification for multidimensional state.
    """
    dims = 2
    inp = [0, 1]
    steps = 8
    constraint_degree = 4
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    # This is a place filler
    constants = [[]] * steps
    extension_factor = 8
    constraint_degree = 4
    def step_fn(f, prev, constants):
      f_n_minus_1 = prev[0]
      f_n = prev[1]
      f_n_plus_1 = f.add(f_n, f_n_minus_1)
      return [f_n, f_n_plus_1]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    proof = mk_proof(comp, params)
    assert verify_proof(comp, params, proof)

  def test_computation_polynomial(self):
    """
    Tests construction of computation polynomial
    """
    dims = 1
    inp = [5]
    steps = 512
    extension_factor = 8
    constants = [[(i**7) ^ 42] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    extension_factor = 8
    constraint_degree = 4
    def step_fn(f, state, constants):
      # 2value**2 + constant
      return [f.add(f.mul(f.exp(state[0], 2), 2), constants[0])]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    assert len(comp_poly_evals) == steps * extension_factor

  def test_constraint_polynomial(self):
    """
    Tests construction of constraint polynomial.
    """
    dims = 1
    inp = [5]
    steps = 512
    extension_factor = 8
    constants = [[(i**7) ^ 42] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    extension_factor = 8
    constraint_degree = 4
    def step_fn(f, state, constants):
      # 2value**2 + constant
      return [f.add(f.mul(f.exp(state[0], 2), 2), constants[0])]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    comp_poly_evals = construct_computation_polynomial(
        comp, params)
    constraint_evals = construct_constraint_polynomial(
        comp, params, comp_poly_evals)
    assert len(constraint_evals) == steps * extension_factor

  def test_compressed_stark(self):
    """Basic compressed stark test"""
    dims = 1
    inp = [3]
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    extension_factor = 8
    constraint_degree = 4
    # Full STARK test
    constants = [[(i**7) ^ 42] for i in range(64)]
    constants = constants * (steps // 64)
    # Factoring out computation
    def step_fn(f, state, constants):
      return [f.add(f.exp(state[0], 3), constants[0])]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
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
    inp = [5]
    steps = 512
    constraint_degree = 4
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    constants = [[(i**7) ^ 42] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    extension_factor = 8

    # Factoring out computation
    def step_fn(f, state, constants):
      # 2value**2 + constant
      return [f.add(f.mul(f.exp(state[0], 2), 2), constants[0])]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
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
    inp = [5]
    steps = 512
    constraint_degree = 4
    constants = [[(i**7) ^ 42] for i in range(64)]
    constants = constants * (steps // 64) 
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    extension_factor = 8

    def step_fn(f, state, constants):
      return [f.add(f.exp(state[0], 3), constants[0])]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    proof = mk_proof(comp, params)
    result = verify_proof(comp, params, proof)
    assert result

  def test_affine_stark(self):
    """
    Basic tests of affine stark generation
    """
    dims = 1
    inp = [5]
    steps = 512
    constraint_degree = 4
    # TODO(rbharath): Why do these constants make sense? Read
    # MiMC paper to see if justification.
    constants = [[0]] * 512
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    extension_factor = 8

    # Factoring out computation
    def step_fn(f, state, constants):
      return [f.add(f.add(f.mul(3, state[0]), 4), constants[0])]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
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
    inp = [5]
    steps = 512
    constraint_degree = 4
    constants = [[i, i] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = PrimeField(modulus)
    ## Factoring out computation
    def step_fn(f, state, constants):
      # c_1*value**2 + c_0
      value = state[0]
      return [f.add(f.mul(constants[1], f.exp(value, 2)), constants[0])]
    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)


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
    inp = [5]
    dims = 1
    steps = 512
    constraint_degree = 4
    constants = [[i, i] for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    extension_factor = 8

    ## Factoring out computation
    def step_fn(f, value, constants):
      # c_1*value**2 + c_0
      return [f.add(f.mul(constants[1], f.exp(value[0], 2)), constants[0])]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    proof = mk_proof(comp, params)
    assert isinstance(proof, list)
    assert len(proof) == 4
    result = verify_proof(comp, params, proof)
    assert result

  def test_varying_quintic_stark(self):
    """
    Basic tests of quintic stark generation
    """
    inp = [5]
    steps = 512
    dims = 1
    constraint_degree = 8
    constants = [[i]*6 for i in range(steps)]
    modulus = 2**256 - 2**32 * 351 + 1
    field = PrimeField(modulus)
    extension_factor = 8

    ## Factoring out computation
    def step_fn(f, state, constants):
      # c_5*value**5 + c_4*value**4 + c_3*value**3 + c_2*value**2 + c_1*value**1 + c_0
      value = state[0]
      return [f.add(f.mul(constants[5], f.exp(value, 5)),
          f.add(f.mul(constants[4], f.exp(value, 4)),
            f.add(f.mul(constants[3], f.exp(value, 3)),
              f.add(f.mul(constants[2], f.exp(value, 2)),
                f.add(f.mul(constants[1], f.exp(value, 1)), constants[0])))))]

    comp = Computation(field, dims, inp, steps, constants, step_fn,
        constraint_degree, extension_factor)
    params = StarkParams(field, comp, modulus, extension_factor)
    proof = mk_proof(comp, params)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    trace, output = get_computational_trace(
        field, inp, steps, constants, step_fn)
    result = verify_proof(comp, params, proof)
    assert result
