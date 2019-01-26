import unittest
from starks.fri import AffineSubspaceFRI
from starks.fri import SmoothSubgroupFRI
from starks.merkle_tree import merkelize
from starks.merkle_tree import merkelize_polynomial_evaluations
from starks.compression import bin_length
from starks.compression import compress_fri
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.fft import NonBinaryFFT
from starks.poly_utils import multivariates_over
from starks.air import Computation
from starks.stark import get_power_cycle 
from starks.stark import STARK
from starks.stark import construct_trace_polynomials
from starks.stark import construct_constraint_polynomials
from starks.stark import construct_remainder_polynomials
from starks.stark import construct_boundary_polynomials
from starks.stark import compute_pseudorandom_linear_combination


class TestFRI(unittest.TestCase):
  """
  Basic tests for FRI implementation. 
  """

  def test_basic_prove(self):
    """Test proof on low degree implementation"""
    degree = 4
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    polysOver = polynomials_over(field).factory
    # 1 + x + 3x^2 + 4 x^3 mod 31
    poly = polysOver([val for val in range(degree)])

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = field(7)**((modulus-1)//8)

    # This is a low degree polynomial so we hit the special
    # case of the handler.
    fri = FRI(field)
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree)
    # The proof is a list of length one, whose first entry is just the evaluations converted to bytes
    assert len(proof[0]) == 8 

  def test_binary_fri_proof(self):
    """Test proof on low degree implementation"""
    # This finite field is of size 2^17
    p = 2
    m = 17
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    field = FiniteField(p, m)
    #field = FiniteField(p, m)
    #x^17 + x^3 + 1 is primitive 
    coefficients = [Zp(0)] * 18
    coefficients[0] = Zp(1)
    coefficients[3] = Zp(1)
    coefficients[17] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    polysOver = polynomials_over(field).factory
    # 1 + x + 3x^2 + 4 x^3 mod 31
    poly = polysOver([val for val in range(degree)])

    # This is a low degree polynomial so we hit the special
    # case of the handler.
    fri = AffineSubspaceFRI(field)
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree)
    # The proof is a list of length one, whose first entry is just the evaluations converted to bytes
    assert len(proof[0]) == 8 

  def test_high_degree_prove(self):
    """Tests proof generation on high degree polynomials"""
    steps = 512 
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    polysOver = polynomials_over(field).factory
    # Some round constants borrowed from MiMC
    poly = polysOver([field((i**7) ^ 42) for i in range(steps)])
    # Root of unity such that x^steps=1
    root_of_unity = field(7)**((modulus-1)//steps)
    # We're trying to prove this is a (steps-1)-degree
    # polnomial
    # degree = (steps-1) + 1 = steps
    fri = FRI(field)
    degree = steps
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree)
    # The proof recurses by dividing maxdeg_plus_1 by 4
    # So 512, 128, 32, 8. (The base case passes over to
    # special handler for degree 16 or less so these are all
    # recursions).
    assert len(proof) == 4
    for i, rec_proof in enumerate(proof):
      if i < 3:
        # Each subproof is [merkle_root, branches] for all but
        # base case.
        assert len(rec_proof) == 2
        assert len(rec_proof[1]) == 40
      else:
        # Here we trigger the base case.
        assert len(rec_proof) == 8

  def test_verify_low_degree_proof(self):
    """Verify a low degree proof"""
    dims = 1
    modulus = 31
    steps = 512
    degree = steps
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    polysOver = polynomials_over(field).factory
    poly = polysOver([field((i**7) ^ 42) for i in range(steps)])
    # Root of unity such that x^steps=1
    root_of_unity = field(7)**((modulus - 1) // steps)
    fri = FRI(field)
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree)

    # TODO(rbharath): Should this be a method?
    fft_solver = NonBinaryFFT(field, root_of_unity)
    evaluations = fft_solver.fft(poly)
    e_mtree = merkelize(evaluations)
    mroot = e_mtree[1]
    verification = fri.verify_proximity_proof(proof, mroot, root_of_unity, degree)
    assert verification

  def test_fri(self):
    """Pure FRI tests"""
    #degree = 4096
    degree = 256
    modulus = 2**256 - 2**32 * 351 + 1
    field= IntegersModP(modulus)
    polysOver = polynomials_over(field).factory
    poly = polysOver([field(val) for val in range(degree)])
    root_of_unity = field(7)**((modulus - 1) // (degree*4))
    #evaluations = fft(poly, modulus, root_of_unity)
    #evaluations = [val[0] for val in evaluations]
    #proof = prove_low_degree(evaluations, root_of_unity, 4096, modulus)
    fri = FRI(field)
    proof = fri.generate_proximity_proof(poly, root_of_unity, degree)
    print("Approx proof length: %d" % bin_length(compress_fri(proof)))
    #assert verify_low_degree_proof(
    #    merkelize(evaluations)[1], root_of_unity, proof, 4096, modulus)

    fft_solver = NonBinaryFFT(field, root_of_unity)
    evaluations = fft_solver.fft(poly)
    e_mtree = merkelize(evaluations)
    mroot = e_mtree[1]
    verification = fri.verify_proximity_proof(proof, mroot, root_of_unity, degree)
    assert verification

    # TODO(rbharath): Make a good test for failure of high degree polynomials
    #fakedata = [
    #    x if pow(3, i, degree) > 400 else 39 for x, i in enumerate(evaluations)
    #]
    #proof2 = prove_low_degree(fakedata, root_of_unity, 4096, modulus)
    #assert not verify_low_degree_proof(
    #    merkelize(fakedata)[1], root_of_unity, proof, 4096, modulus)

    #try:
    #  assert verify_low_degree_proof(
    #      merkelize(evaluations)[1], root_of_unity, proof, 2048, modulus)
    #  raise Exception("Fake data passed FRI")
    #except:
    #  pass

  def test_varying_quadratic_fri(self):
    """
    Basic tests of FRI generation for quadratic stark with varying coefficients
    """
    width = 2
    steps = 2
    modulus = 2**256 - 2**32 * 351 + 1
    extension_factor = 8 
    field = IntegersModP(modulus)
    inp = [field(2), field(5)]

    precision = steps * extension_factor
    G1 = field(7)**((modulus - 1) // steps)
    G2 = field(7)**((modulus - 1) // precision)
    xs = get_power_cycle(G2, field)
    last_step_position = xs[(steps - 1) * extension_factor]
    ### Factoring out computation
    polysOver = multivariates_over(field, width).factory
    X_1 = polysOver({(1,0): field(1)})
    X_2 = polysOver({(0,1): field(1)})
    step_polys = [X_1, X_1 + X_2**2] 
    comp = Computation(field, width, inp, steps, step_polys,
        extension_factor)
    witness = comp.generate_witness()
    boundary = comp.generate_boundary_constraints()
    stark = STARK(field, steps, modulus, extension_factor, width, step_polys)


    trace_polys = construct_trace_polynomials(witness, field, G1)
    constraint_polys = construct_constraint_polynomials(step_polys,
        trace_polys, field, G1, width)
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
    l_root = l_mtree[1]

    fri = FRI(field)
    proof = fri.generate_proximity_proof(l_poly, G2, steps*comp.get_degree(), exclude_multiples_of=extension_factor)

    fft_solver = NonBinaryFFT(field, G2)
    evaluations = fft_solver.fft(l_poly)
    e_mtree = merkelize(evaluations)
    mroot = e_mtree[1]
    verification = fri.verify_proximity_proof(proof, mroot, G2, steps*comp.get_degree(), exclude_multiples_of=extension_factor)
    assert verification
