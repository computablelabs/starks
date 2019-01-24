import unittest
from starks.merkle_tree import merkelize
from starks.merkle_tree import merkelize_polynomial_evaluations
from starks.air import Computation
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
from starks.poly_utils import multivariates_over


class TestStark(unittest.TestCase):
  """
  Basic tests for Stark construction implementation. 
  """


  def test_stark_init(self):
    """Test generation of stark parameters."""
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    extension_factor = 8
    # Only tests that constructor works implicitly
    params = STARK(field, steps, modulus, extension_factor, width=1, step_polys=[])

  def test_binary_stark_init(self):
    """Test generation of stark parameters."""
    steps = 512
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    extension_factor = 8
    # Only tests that constructor works implicitly
    params = STARK(field, steps, modulus, extension_factor, width=1, step_polys=[])


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

  def test_higher_dim_proof_verification(self):
    """
    Tests proof generation and verification for multidimensional state.
    """
    width = 2
    steps = 32 
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(0), field(1)]
    extension_factor = 8
    [X_1, X_2] = generate_Xi_s(field, width)
    step_polys = [X_2, X_1 + X_2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()

    stark = STARK(field, steps, modulus, extension_factor, width, step_polys)
    proof = stark.mk_proof(witness, boundary)
    assert stark.verify_proof(proof, witness, boundary)

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

    stark = STARK(field, steps, modulus, extension_factor, width, step_polys)
    proof = stark.mk_proof(witness, boundary)
    result = stark.verify_proof(proof, witness, boundary)
    assert result

  def test_affine_stark(self):
    """
    Basic tests of affine stark generation
    """
    width = 2
    steps = 32 
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    inp = [field(2), field(5)]
    extension_factor = 8

    ## Factoring out computation
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_1, X_1 + 3*X_2] 

    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()

    stark = STARK(field, steps, modulus, extension_factor, width, step_polys)
    proof = stark.mk_proof(witness, boundary)
    assert isinstance(proof, list)
    assert len(proof) == 4
    result = stark.verify_proof(proof, witness, boundary)
    assert result

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

    stark = STARK(field, steps, modulus, extension_factor, width, step_polys)
    proof = stark.mk_proof(witness, boundary)
    assert isinstance(proof, list)
    assert len(proof) == 4
    (m_root, l_root, branches, fri_proof) = proof
    result = stark.verify_proof(proof, witness, boundary)
    assert result
