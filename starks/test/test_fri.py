import math
import unittest
from starks.fri import AffineSubspaceFRI
#from starks.fri import SmoothSubgroupFRI
from starks.merkle_tree import merkelize
from starks.merkle_tree import merkelize_polynomial_evaluations
from starks.compression import bin_length
from starks.compression import compress_fri
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.fft import NonBinaryFFT
from starks.poly_utils import multivariates_over
from starks.air import AIR
from starks.utils import get_power_cycle
#from starks.stark import construct_trace_polynomials
#from starks.stark import construct_constraint_polynomials
#from starks.stark import construct_remainder_polynomials
#from starks.stark import construct_boundary_polynomials
#from starks.stark import compute_pseudorandom_linear_combination
from starks.finitefield import FiniteField
from starks.reedsolomon import AffineSpace
from starks.fri import select_subspace
from starks.fri import select_random
from starks.fri import eval_on_subspace
from starks.poly_utils import construct_affine_vanishing_polynomial
from starks.poly_utils import lagrange_interp


class TestFRI(unittest.TestCase):
  """
  Basic tests for FRI implementation.
  """

  #def test_basic_prove(self):
  #  """Test proof on low degree implementation"""
  #  degree = 4
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  field = IntegersModP(modulus)
  #  polysOver = polynomials_over(field).factory
  #  # 1 + x + 3x^2 + 4 x^3 mod 31
  #  poly = polysOver([val for val in range(degree)])

  #  # A root of unity is a number such that z^n = 1
  #  # This provides us a 6-th root of unity (z^6 = 1)
  #  root_of_unity = field(7)**((modulus-1)//8)

  #  # This is a low degree polynomial so we hit the special
  #  # case of the handler.
  #  fri = SmoothSubgroupFRI(field)
  #  proof = fri.generate_proximity_proof(poly, root_of_unity, degree)
  #  # The proof is a list of length one, whose first entry is just the evaluations converted to bytes
  #  assert len(proof[0]) == 8

  def test_binary_fri_commit_query(self):
    """Test commit and query in Fri"""
    
    # testing commit part
    p = 2
    m = 10
    Zp = IntegersModP(p)
    basePolys = polynomials_over(Zp)
    g = basePolys([0, 1])
    field = FiniteField(p, m)
    t = 2
    S = AffineSpace(Zp, [g**k for k in range(t)])
    rho = 1   
    fri = AffineSubspaceFRI(field, S, rho)
    fri_poly = basePolys([Zp(0), Zp(1)])
    vals = fri.commit(fri_poly)
    print("vals")
    print(vals)


    #testing query part
    R = -math.log(rho, 2)
    eta = 2
    k_0 = math.log(len(S), 2)
    r = int(math.floor((k_0 - R)/eta))
    L_i = S
    L_i_0 = select_subspace(L_i, eta)
    f = []
    f.append(fri_poly)
    for i in range(r-1):
        x = str(R + eta + r + k_0 + i)
        x_i = select_random(field, x)
        q_i = construct_affine_vanishing_polynomial(field, L_i_0)
        L_i_plus_1 = eval_on_subspace(q_i, L_i)
        f_i_plus_1 = {}
        for y in L_i_plus_1:
            coset = [y + l_i for l_i in L_i_0]
            coset_eval = [f[-1](x_c) for x_c in coset]
            interp = lagrange_interp(self.field, coset, coset_eval)
            f_i_plus_1 = interp(x_i)
        f.append(lambda z: f_i_plus_1[z])
        L_i = L_i_plus_1
        L_i_0 = select_subspace(L_i, eta)

    result_verify = fri.query(f, vals)
    assert result_verify == True

#    proof = fri.generate_proximity_proof(poly, S)
#    # The proof is a list of length one, whose first entry is just the evaluations converted to bytes
#    assert len(proof[0]) == 8

  #def test_high_degree_prove(self):
  #  """Tests proof generation on high degree polynomials"""
  #  steps = 512
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  field = IntegersModP(modulus)
  #  polysOver = polynomials_over(field).factory
  #  # Some round constants borrowed from MiMC
  #  poly = polysOver([field((i**7) ^ 42) for i in range(steps)])
  #  # Root of unity such that x^steps=1
  #  root_of_unity = field(7)**((modulus-1)//steps)
  #  # We're trying to prove this is a (steps-1)-degree
  #  # polnomial
  #  # degree = (steps-1) + 1 = steps
  #  fri = SmoothSubgroupFRI(field)
  #  degree = steps
  #  proof = fri.generate_proximity_proof(poly, root_of_unity, degree)
  #  # The proof recurses by dividing maxdeg_plus_1 by 4
  #  # So 512, 128, 32, 8. (The base case passes over to
  #  # special handler for degree 16 or less so these are all
  #  # recursions).
  #  assert len(proof) == 4
  #  for i, rec_proof in enumerate(proof):
  #    if i < 3:
  #      # Each subproof is [merkle_root, branches] for all but
  #      # base case.
  #      assert len(rec_proof) == 2
  #      assert len(rec_proof[1]) == 40
  #    else:
  #      # Here we trigger the base case.
  #      assert len(rec_proof) == 8

  #def test_verify_low_degree_proof(self):
  #  """Verify a low degree proof"""
  #  dims = 1
  #  modulus = 31
  #  steps = 512
  #  degree = steps
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  field = IntegersModP(modulus)
  #  polysOver = polynomials_over(field).factory
  #  poly = polysOver([field((i**7) ^ 42) for i in range(steps)])
  #  # Root of unity such that x^steps=1
  #  root_of_unity = field(7)**((modulus - 1) // steps)
  #  fri = SmoothSubgroupFRI(field)
  #  proof = fri.generate_proximity_proof(poly, root_of_unity, degree)

  #  # TODO(rbharath): Should this be a method?
  #  fft_solver = NonBinaryFFT(field, root_of_unity)
  #  evaluations = fft_solver.fft(poly)
  #  e_mtree = merkelize(evaluations)
  #  mroot = e_mtree[1]
  #  verification = fri.verify_proximity_proof(proof, mroot, root_of_unity, degree)
  #  assert verification

  #def test_fri(self):
  #  """Pure FRI tests"""
  #  #degree = 4096
  #  degree = 256
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  field= IntegersModP(modulus)
  #  polysOver = polynomials_over(field).factory
  #  poly = polysOver([field(val) for val in range(degree)])
  #  root_of_unity = field(7)**((modulus - 1) // (degree*4))
  #  #evaluations = fft(poly, modulus, root_of_unity)
  #  #evaluations = [val[0] for val in evaluations]
  #  #proof = prove_low_degree(evaluations, root_of_unity, 4096, modulus)
  #  fri = SmoothSubgroupFRI(field)
  #  proof = fri.generate_proximity_proof(poly, root_of_unity, degree)
  #  print("Approx proof length: %d" % bin_length(compress_fri(proof)))
  #  #assert verify_low_degree_proof(
  #  #    merkelize(evaluations)[1], root_of_unity, proof, 4096, modulus)

  #  fft_solver = NonBinaryFFT(field, root_of_unity)
  #  evaluations = fft_solver.fft(poly)
  #  e_mtree = merkelize(evaluations)
  #  mroot = e_mtree[1]
  #  verification = fri.verify_proximity_proof(proof, mroot, root_of_unity, degree)
  #  assert verification

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

  #def test_varying_quadratic_fri(self):
  #  """
  #  Basic tests of FRI generation for quadratic stark with varying coefficients
  #  """
  #  width = 2
  #  steps = 3
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  extension_factor = 8
  #  field = IntegersModP(modulus)
  #  inp = [field(2), field(5)]

  #  precision = steps * extension_factor
  #  G1 = field(7)**((modulus - 1) // steps)
  #  G2 = field(7)**((modulus - 1) // precision)
  #  xs = get_power_cycle(G2, field)
  #  last_step_position = xs[(steps - 1) * extension_factor]
  #  ### Factoring out computation
  #  polysOver = multivariates_over(field, width).factory
  #  X_1 = polysOver({(1,0): field(1)})
  #  X_2 = polysOver({(0,1): field(1)})
  #  step_polys = [X_1, X_1 + X_2**2]
  #  comp = AIR(field, width, inp, steps, step_polys,
  #      extension_factor)
  #  witness = comp.generate_witness()
  #  boundary = comp.generate_boundary_constraints()
  #  stark = STARK(field, steps, modulus, extension_factor, width, step_polys)


  #  trace_polys = construct_trace_polynomials(witness, field, G1)
  #  constraint_polys = construct_constraint_polynomials(step_polys,
  #      trace_polys, field, G1, width)
  #  remainder_polys = construct_remainder_polynomials(constraint_polys, field,
  #      steps, last_step_position)
  #  boundary_polys = construct_boundary_polynomials(
  #      trace_polys, witness, boundary, field, last_step_position, width)

  #  fft_solver = NonBinaryFFT(field, G2)
  #  polys = trace_polys + remainder_polys + boundary_polys
  #  poly_evals = []
  #  for poly in polys:
  #    poly_eval = fft_solver.fft(poly)
  #    poly_evals.append(poly_eval)
  #  mtree = merkelize_polynomial_evaluations(width, poly_evals)

  #  l_poly = compute_pseudorandom_linear_combination(mtree[1],
  #      trace_polys, remainder_polys, boundary_polys, field, G2, precision,
  #      steps, width)
  #  l_evaluations = fft_solver.fft(l_poly)
  #  l_mtree = merkelize(l_evaluations)
  #  l_root = l_mtree[1]

  #  fri = FRI(field)
  #  proof = fri.generate_proximity_proof(l_poly, G2, steps*comp.get_degree(), exclude_multiples_of=extension_factor)

  #  fft_solver = NonBinaryFFT(field, G2)
  #  evaluations = fft_solver.fft(l_poly)
  #  e_mtree = merkelize(evaluations)
  #  mroot = e_mtree[1]
  #  verification = fri.verify_proximity_proof(proof, mroot, G2, steps*comp.get_degree(), exclude_multiples_of=extension_factor)
  #  assert verification
