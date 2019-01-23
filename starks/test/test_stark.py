import unittest
import time
from starks.utils import mimc
from starks.merkle_tree import merkelize
from starks.merkle_tree import verify_branch
from starks.merkle_tree import mk_branch
from starks.merkle_tree import blake
from starks.merkle_tree import merkelize_polynomial_evaluations
from starks.merkle_tree import unpack_merkle_leaf
from starks.air import Computation
from starks.air import get_computational_trace
from starks.fft import NonBinaryFFT
# TODO(rbharath): These need to be swapped out for correct imports
from starks.utils import generate_Xi_s
from starks.utils import get_pseudorandom_indices
from starks.stark import get_power_cycle 
from starks.stark import construct_trace_polynomials
from starks.stark import construct_constraint_polynomials
from starks.stark import construct_remainder_polynomials
from starks.stark import construct_boundary_polynomials
from starks.stark import get_pseudorandom_ks
from starks.stark import compute_pseudorandom_linear_combination_1d
from starks.stark import compute_pseudorandom_linear_combination
from starks.stark import STARK 
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

  def test_trace_polynomials(self):
    """
    Tests construction of computation polynomial
    """
    width = 2
    steps = 128 
    extension_factor = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(2), field(5)]
    [X_1, X_2] = generate_Xi_s(field, width)
    step_polys = [X_2, X_1 + 2*X_2**2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    witness = comp.generate_witness()
    params = STARK(field, steps, modulus, extension_factor, width, step_polys)
    trace_polys = construct_trace_polynomials(witness, params.field, params.G1)
    assert len(trace_polys) == width
    xs = get_power_cycle(params.G1, params.field) 
    # Check that the trace polynomial reconstitutes the witness
    for dim in range(width):
      witness_dim = witness[dim]
      trace_poly = trace_polys[dim]
      for ind, x in enumerate(xs):
        assert witness_dim[ind] == trace_poly(x)

  def test_constraint_polynomials(self):
    """
    Tests construction of constraint polynomial.
    """
    width = 2
    steps = 4 
    extension_factor = 2
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    G1 = field(7)**((modulus - 1) // steps)
    inp = [field(2), field(5)]
    [X_1, X_2] = generate_Xi_s(field, width)
    step_polys = [X_2, X_1 + 2*X_2**2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    witness = comp.generate_witness()

    trace_polys = construct_trace_polynomials(witness, field, G1)
    constraint_polys = construct_constraint_polynomials(step_polys, trace_polys, field, G1, width)
    assert len(constraint_polys) == width
    # TODO(rbharath): Add more meaningful test here.

  def test_pseudorandom_combo_merkle_root(self):
    """
    Tests that merkle spot checks are constructed correctly.
    """
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    width = 3
    # state = [c_0, c_1, value]
    inp = [field(2), field(2), field(5)]
    steps = 4
    spot_check_security_factor = 80
    extension_factor = 8
    precision = steps * extension_factor
    G1 = field(7)**((modulus - 1) // steps)
    G2 = field(7)**((modulus - 1) // precision)
    xs = get_power_cycle(G2, field)
    last_step_position = xs[(steps - 1) * extension_factor]

    [X_1, X_2, X_3] = generate_Xi_s(field, width)
    step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
    comp = Computation(field, width, inp, steps, step_polys, extension_factor)
    #params = STARK(field, steps, modulus, extension_factor, width, step_polys)

    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()

    trace_polys = construct_trace_polynomials(witness, field, G1)
    constraint_polys = construct_constraint_polynomials(step_polys, trace_polys, field, G1, width)
    remainder_polys = construct_remainder_polynomials(constraint_polys, field,
        steps, last_step_position)
    boundary_polys = construct_boundary_polynomials(
        trace_polys, witness, boundary, field, last_step_position, width)

    fft_solver = NonBinaryFFT(field, G2)
    polys = trace_polys + remainder_polys + boundary_polys
    poly_evals = []
    for poly in polys:
      poly_eval = fft_solver.fft(poly)
      poly_evals.append(poly_eval)
    mtree = merkelize_polynomial_evaluations(width, poly_evals)
    l_poly = compute_pseudorandom_linear_combination(mtree[1],
        trace_polys, remainder_polys, boundary_polys, field, G2, precision,
        steps, width)
    l_evaluations = fft_solver.fft(l_poly)
    l_mtree = merkelize(l_evaluations)

    # TODO(rbharath): Move this to anouther test
    ## Compute merkle spot checks
    #branches = stark.compute_merkle_spot_checks(mtree, l_mtree,
    #    samples=spot_check_security_factor)

    ## Each spot check returns 3 branches, m[pos], m[pos+extension_factor], l[pos] 
    #assert len(branches) == 3*spot_check_security_factor

  def test_stark_init(self):
    """Test generation of stark parameters."""
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    extension_factor = 8
    # Only tests that constructor works implicitly
    params = STARK(field, steps, modulus, extension_factor, width=1, step_polys=[])

  def test_get_pseudorandom_indices(self):
    """
    Tests that pseudorandom elements are computed correctly.
    """
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(2), field(2), field(5)]
    width = 3
    steps = 4
    spot_check_security_factor = 80
    polysOver = multivariates_over(field, width).factory
    extension_factor = 8
    ## Factoring out computation
    [X_1, X_2, X_3] = generate_Xi_s(field, width)
    step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
    precision = steps * extension_factor
    G1 = field(7)**((modulus - 1) // steps)
    G2 = field(7)**((modulus - 1) // precision)
    xs = get_power_cycle(G2, field)
    last_step_position = xs[(steps - 1) * extension_factor]

    comp = Computation(field, width, inp, steps, step_polys, extension_factor)
    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()

    trace_polys = construct_trace_polynomials(witness, field, G1)
    constraint_polys = construct_constraint_polynomials(step_polys, trace_polys, field, G1, width)
    remainder_polys = construct_remainder_polynomials(constraint_polys, field,
        steps, last_step_position)
    boundary_polys = construct_boundary_polynomials(
        trace_polys, witness, boundary, field, last_step_position, width)

    fft_solver = NonBinaryFFT(field, G2)
    polys = trace_polys + remainder_polys + boundary_polys
    poly_evals = []
    for poly in polys:
      poly_eval = fft_solver.fft(poly)
      poly_evals.append(poly_eval)
    mtree = merkelize_polynomial_evaluations(width, poly_evals)
    l_poly = compute_pseudorandom_linear_combination(mtree[1],
        trace_polys, remainder_polys, boundary_polys, field, G2, precision,
        steps, width)
    l_evaluations = fft_solver.fft(l_poly)
    l_mtree = merkelize(l_evaluations)

    indices = get_pseudorandom_indices(l_mtree[1],
        precision, count=spot_check_security_factor,
        exclude_multiples_of=extension_factor)
    assert len(indices) == spot_check_security_factor

  def test_verify_pseudorandom_combination(self):
    """
    Tests compute pseudorandom linear combination for 1-dimension
    """
    width = 3
    steps = 8
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8
    field = IntegersModP(modulus)
    # state = [c_0, c_1, value]
    inp = [field(2), field(2), field(5)]
    ## Factoring out computation
    [X_1, X_2, X_3] = generate_Xi_s(field, width)
    step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
    precision = steps * extension_factor
    G1 = field(7)**((modulus - 1) // steps)
    G2 = field(7)**((modulus - 1) // precision)
    xs = get_power_cycle(G2, field)
    last_step_position = xs[(steps - 1) * extension_factor]

    comp = Computation(field, width, inp, steps, step_polys, extension_factor)
    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()

    trace_polys = construct_trace_polynomials(witness, field, G1)
    constraint_polys = construct_constraint_polynomials(step_polys, trace_polys, field, G1, width)
    remainder_polys = construct_remainder_polynomials(constraint_polys, field,
        steps, last_step_position)
    boundary_polys = construct_boundary_polynomials(
        trace_polys, witness, boundary, field, last_step_position, width)

    fft_solver = NonBinaryFFT(field, G2)
    trace_evals = []
    for poly in trace_polys:
      trace_eval = fft_solver.fft(poly)
      trace_evals.append(trace_eval)
    remainder_evals = []
    for poly in remainder_polys:
      remainder_eval = fft_solver.fft(poly)
      remainder_evals.append(remainder_eval)
    boundary_evals = []
    for poly in boundary_polys:
      boundary_eval = fft_solver.fft(poly)
      boundary_evals.append(boundary_eval)
    poly_evals = trace_evals + remainder_evals + boundary_evals
    mtree = merkelize_polynomial_evaluations(width, poly_evals)
    
    entropy = mtree[1]
    l_polys_per_dim = compute_pseudorandom_linear_combination_1d(
        entropy, trace_polys, remainder_polys, boundary_polys, G2, steps, precision)

    for dim in range(width):
      for i, pos in enumerate(range(precision)):
        p_of_x = trace_evals[dim][pos]
        next_pos = (pos + extension_factor) % precision
        p_of_g1x = trace_evals[dim][next_pos]
        d_of_x = remainder_evals[dim][pos]
        b_of_x = boundary_evals[dim][pos]

        x = G2**pos
        x_to_the_steps = x**steps

        # Leaf node from l[pos]
        k1, k2, k3, k4 = get_pseudorandom_ks(mtree[1], 4)
        for dim in range(width):
          l_of_x_dim = l_polys_per_dim[dim](x)
          assert (l_of_x_dim - d_of_x - k1 * p_of_x - k2 * p_of_x * x_to_the_steps - k3 * b_of_x - k4 * b_of_x * x_to_the_steps) == 0

  #def test_1d_end_to_end(self):
  #  """
  #  Tests stark end-to-end for 1d example 
  #  """
  #  width = 3
  #  steps = 512
  #  spot_check_security_factor = 80
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  extension_factor = 8
  #  field = IntegersModP(modulus)
  #  inp = [field(2), field(2), field(5)]
  #  ## Factoring out computation
  #  [X_1, X_2, X_3] = generate_Xi_s(field, width)
  #  step_polys = [X_1, X_2, X_1 + X_2*X_3**2] 
  #  comp = Computation(field, width, inp, steps, step_polys,
  #      extension_factor)
  #  params = StarkParams(field, steps, modulus, extension_factor)

  #  p_evaluations = construct_computation_polynomial(
  #      comp, params)
  #  c_of_p_evaluations = construct_constraint_polynomial(
  #      comp, params, p_evaluations)
  #  d_evaluations = construct_remainder_polynomial(
  #      comp, params, c_of_p_evaluations)
  #  b_evaluations = construct_boundary_polynomial(
  #      comp, params, p_evaluations)

  #  polys = [p_evaluations, d_evaluations, b_evaluations]
  #  mtree = merkelize_polynomials(width, polys)
  #  l_evaluations = compute_pseudorandom_linear_combination(
  #      comp, params, mtree, polys)
  #  l_mtree = merkelize(l_evaluations)
  #  m_root = mtree[1]
  #  l_root = l_mtree[1]
  #  fri_proof = prove_low_degree(
  #        l_evaluations,
  #        params.G2,
  #        # Manually filling in degree here
  #        steps * 3,
  #        modulus,
  #        exclude_multiples_of=extension_factor)

  #  assert verify_low_degree_proof(
  #      l_root,
  #      params.G2,
  #      fri_proof,
  #      # Manually filling in degree here
  #      steps * 3,
  #      modulus,
  #      exclude_multiples_of=extension_factor)

  #  # TODO(rbharath): Factor this into function?
  #  byte_list = [b'0x01', b'0x02', b'0x03', b'0x04']
  #  ks = [int.from_bytes(blake(mtree[1] + byte_list[ind]), 'big') for ind in range(4)]
  #  k1, k2, k3, k4 = ks

  #  branches = compute_merkle_spot_checks(mtree, l_mtree, comp, params)
  #  positions = get_pseudorandom_indices(l_mtree[1],
  #      params.precision, count=spot_check_security_factor,
  #      exclude_multiples_of=params.extension_factor)

  #  for i, pos in enumerate(positions):
  #    branch1 = mk_branch(mtree, pos)
  #    mbranch1 = verify_branch(m_root, pos, branch1)
  #    unpacked_leaf1 = unpack_merkle_leaf(mbranch1, comp.width, 3)
  #    next_pos = (pos + params.extension_factor) % params.precision
  #    branch2 = mk_branch(mtree, next_pos) 
  #    mbranch2 = verify_branch(m_root, next_pos, branch2)
  #    unpacked_leaf2 = unpack_merkle_leaf(mbranch2, comp.width, 3)
  #    # Leaf node from l[pos]
  #    l_of_x = verify_branch(l_root, pos, branches[i * 3 + 2],
  #        output_as_int=True)

  #    # Check that p_of_x was recovered correctly
  #    p_of_x = p_evaluations[pos]
  #    p_of_x_recovered = [field(p_of_x_dim) for p_of_x_dim in unpacked_leaf1[:comp.width]]
  #    assert len(p_of_x) == len(p_of_x_recovered)
  #    for dim in range(comp.width):
  #      assert p_of_x[dim] == p_of_x_recovered[dim]

  #    # Check that p_of_g1x was recovered correctly
  #    p_of_g1x = p_evaluations[next_pos]
  #    p_of_g1x_recovered = [field(p_of_g1x_dim) for p_of_g1x_dim in unpacked_leaf2[:comp.width]]
  #    assert len(p_of_g1x) == len(p_of_g1x_recovered)
  #    for dim in range(comp.width):
  #      assert p_of_g1x[dim] == p_of_g1x_recovered[dim]

  #    d_of_x = d_evaluations[pos]
  #    d_of_x_recovered = [field(d_of_x_dim) for d_of_x_dim in unpacked_leaf1[comp.width:2*comp.width]]
  #    assert len(d_of_x) == len(d_of_x_recovered)
  #    for dim in range(comp.width):
  #      assert d_of_x[dim] == d_of_x_recovered[dim]

  #    b_of_x = b_evaluations[pos]
  #    b_of_x_recovered = [field(b_of_x_dim) for b_of_x_dim in unpacked_leaf1[2*comp.width:]]
  #    assert len(b_of_x) == len(b_of_x_recovered)
  #    for dim in range(comp.width):
  #      assert b_of_x[dim] == b_of_x_recovered[dim]

  #    x = params.G2**pos
  #    x_to_the_steps = x**comp.steps
  #    zvalue = (x**comp.steps - 1)/(x - params.last_step_position)
  #    f_of_p_of_x = [comp.step_polys[i](p_of_x) for i in range(width)]
  #    f_of_p_of_x_recovered = [comp.step_polys[i](p_of_x) for i in range(width)]
  #    assert f_of_p_of_x == f_of_p_of_x_recovered
  #    assert (p_of_g1x[0] - f_of_p_of_x[0] - zvalue * d_of_x[0]) == 0
  #    assert (p_of_g1x_recovered[0] - f_of_p_of_x_recovered[0] - zvalue * d_of_x_recovered[0]) == 0

  #    # TODO(rbharath): How do I do this with multidimensional polynomial?
  #   # zeropoly2 = polysOver([-1, 1]) * polysOver([-params.last_step_position, 1])
  #   # for dim in range(comp.width):
  #   #   interpolant_dim = lagrange_interp_2(modulus, [1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
  #   #   assert (p_of_x[dim] - b_of_x[dim] * zeropoly2(x) - interpolant_dim(x)) == 0

  # TODO(rbharath): This is broken!! Need to fix in future PR
  #def test_higher_dim_proof_verification(self):
  #  """
  #  Tests proof generation and verification for multidimensional state.
  #  """
  #  width = 2
  #  steps = 8
  #  constraint_degree = 4
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  field = IntegersModP(modulus)
  #  inp = [field(0), field(1)]
  #  extension_factor = 8
  #  [X_1, X_2] = generate_Xi_s(field, width)
  #  step_polys = [X_2, X_1 + X_2] 
  #  comp = Computation(field, width, inp, steps, step_polys,
  #      extension_factor)
  #  params = StarkParams(field, steps, modulus, extension_factor)
  #  proof = mk_proof(comp, params)
  #  assert verify_proof(comp, params, proof)

  def test_quadratic_stark(self):
    """
    Basic tests of quadratic stark generation
    """
    width = 2
    steps = 8
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(2), field(5)]
    extension_factor = 8

    # Factoring out computation
    # Polys are: [X_1, X_1 + X_2**2]
    [X_1, X_2] = generate_Xi_s(field, width)
    step_polys = [X_1, X_1 + X_2**3] 

    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    stark = STARK(field, steps, modulus, extension_factor, width, step_polys)

    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()
    proof = stark.mk_proof(witness, boundary)
    assert isinstance(proof, list)
    assert len(proof) == 4
    result = stark.verify_proof(proof, witness, boundary)
    assert result


  def test_mimc_stark_verification(self):
    """
    Basic tests of MiMC stark verification.
    """
    width = 2
    steps = 8
    constraint_degree = 4
    # TODO(rbharath): Should be able to encode these constants as boundary
    # conditions on the tape for the STARK
    #constants = [[(i**7) ^ 42] for i in range(64)]
    #constants = constants * (steps // 64) 
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(2), field(5)]
    extension_factor = 8

    ## Factoring out computation
    [X_1, X_2] = generate_Xi_s(field, width)
    step_polys = [X_1, X_1 + X_2**3] 

    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()

    params = STARK(field, steps, modulus, extension_factor, width, step_polys)
    proof = mk_proof(witness, boundary, params)
    result = verify_proof(proof, witness, boundary, params)
    assert result

  # TODO(rbharath): This is broken!! Need to fix in future PR
  #def test_affine_stark(self):
  #  """
  #  Basic tests of affine stark generation
  #  """
  #  width = 2
  #  steps = 512
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  field = IntegersModP(modulus)
  #  inp = [field(2), field(5)]
  #  extension_factor = 8

  #  ## Factoring out computation
  #  polysOver = multivariates_over(field, width).factory
  #  X_1 = polysOver({(1,0): field(1)})
  #  X_2 = polysOver({(0,1): field(1)})
  #  step_polys = [X_1, X_1 + 3*X_2] 

  #  comp = Computation(field, width, inp, steps, step_polys,
  #      extension_factor)
  #  params = StarkParams(field, steps, modulus, extension_factor)
  #  proof = mk_proof(comp, params)
  #  assert isinstance(proof, list)
  #  assert len(proof) == 4
  #  result = verify_proof(comp, params, proof)
  #  assert result

  # TODO(rbharath): This is broken!! Need to fix in future PR
  #def test_varying_quadratic_fri(self):
  #  """
  #  Basic tests of FRI generation for quadratic stark with varying coefficients
  #  """
  #  width = 2
  #  steps = 512
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  extension_factor = 8
  #  field = IntegersModP(modulus)
  #  inp = [field(2), field(5)]
  #  ### Factoring out computation
  #  polysOver = multivariates_over(field, width).factory
  #  X_1 = polysOver({(1,0): field(1)})
  #  X_2 = polysOver({(0,1): field(1)})
  #  step_polys = [X_1, X_1 + X_2**2] 
  #  comp = Computation(field, width, inp, steps, step_polys,
  #      extension_factor)
  #  params = StarkParams(field, steps, modulus, extension_factor)


  #  p_evaluations = construct_computation_polynomial(
  #      comp, params)
  #  c_of_p_evaluations = construct_constraint_polynomial(
  #      comp, params, p_evaluations)
  #  d_evaluations = construct_remainder_polynomial(
  #      comp, params, c_of_p_evaluations)
  #  b_evaluations = construct_boundary_polynomial(
  #      comp, params, p_evaluations)

  #  polys = [p_evaluations, d_evaluations, b_evaluations]
  #  mtree = merkelize_polynomials(width, polys) 
  #  l_evaluations = compute_pseudorandom_linear_combination(
  #      comp, params, mtree, polys)
  #  l_mtree = merkelize(l_evaluations)
  #  l_root = l_mtree[1]
  #  fri_proof = prove_low_degree(
  #        l_evaluations,
  #        params.G2,
  #        steps * comp.get_degree(),
  #        modulus,
  #        exclude_multiples_of=extension_factor)

  #  assert verify_low_degree_proof(
  #      l_root,
  #      params.G2,
  #      fri_proof,
  #      steps * comp.get_degree(),
  #      modulus,
  #      exclude_multiples_of=extension_factor)

  def test_varying_quintic_stark(self):
    """
    Basic tests of quintic stark generation
    """
    steps = 8 
    width = 6
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(1), field(2), field(3), field(4), field(5), field(6)]
    extension_factor = 8

    ### Factoring out computation
    [X_1, X_2, X_3, X_4, X_5, X_6] = generate_Xi_s(field, width)
    step_polys = [X_1, X_2, X_3, X_4, X_5, X_1*X_2*X_3*X_4*X_5*X_6]

    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()

    params = STARK(field, steps, modulus, extension_factor, width, step_polys)
    proof = mk_proof(witness, boundary, params)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    trace, output = get_computational_trace(
        inp, steps, width, step_polys)
    result = verify_proof(proof, witness, boundary, params)
    assert result
