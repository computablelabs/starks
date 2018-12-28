import time
from starks.merkle_tree import blake
from starks.merkle_tree import verify_branch
from starks.merkle_tree import mk_branch
from starks.merkle_tree import merkelize
from starks.merkle_tree import merkelize_polynomials
from starks.merkle_tree import unpack_merkle_leaf
from starks.polynomial import polynomials_over
from starks.poly_utils import lagrange_interp_2 
from starks.numbertype import FieldElement
from starks.fft import fft
from starks.fri import prove_low_degree, verify_low_degree_proof
from starks.utils import get_power_cycle, get_pseudorandom_indices, is_a_power_of_2
from starks.air import Computation
from starks.poly_utils import multi_inv

# Number of branches used for Merkle-tree check


class StarkParams(object):
  """Holds the cryptographic parameters needed for STARK"""
  def __init__(self, field, steps: int, modulus: FieldElement,
      extension_factor: int, spot_check_security_factor: int =80):
    """
    Parameters
    ----------
    field: Field
      The Field in which computation is permored
    steps: int 
      The number of steps in Computation 
    modulus: Int
      A prime p that defines finite field Z/p
    extension_factor: Int
      A power of two which is the degree to which the trace is expanded
      when  constructing polynomials. For example, a trace of length 512
      with an extension_factor of 8 would construct polynomial evaluations
      on a 4096 elements.
    """
    self.field = field
    self.modulus = modulus
    self.extension_factor = extension_factor
    self.precision = steps * extension_factor
    self.spot_check_security_factor = spot_check_security_factor

    # Root of unity such that x^precision=1
    #self.G2 = self.field.exp(7, (modulus - 1) // self.precision)
    self.G2 = field(7)**((modulus - 1) // self.precision)

    # Root of unity such that x^steps=1
    #self.G1 = self.field.exp(self.G2, extension_factor)
    self.G1 = self.G2**extension_factor

    ## Powers of the higher-order root of unity
    self.xs = get_power_cycle(self.G2, modulus)
    self.last_step_position = self.xs[(steps - 1) * extension_factor]

def construct_constants_polynomials(comp, params):
  """Transforms constants into polynomials
  
  TODO(rbharath): Refactoring the constants list to be a list of step-wise
  constants might be better.

  The constants needed for a computation are encoded as a list of lists.
  Each of these lists is of length comp.steps. These encapsulate the
  external state for the computation at each step. This function
  interpolates these constants lists into polynomials. These are used by
  the prover and the verifier to respectively prove and verify facts about
  computation.
  """
  constants_polynomials = []
  constants_extensions = []
  # This is safe since steps > 0
  deg = len(comp.constants[0])
  for d in range(deg):
    # The extra wrapping is some plumbing since the fft expects a sequence
    # of states, where a state is a list.
    deg_constants = [[step_constants[d]] for step_constants in comp.constants]
    # Constants are a 1-d sequence
    constants_mini_polynomial = fft(
        deg_constants, params.modulus, params.G1,
        inv=True, dims=1)
    constants_mini_extension = fft(constants_mini_polynomial,
        params.modulus, params.G2, dims=1)
    assert len(constants_mini_extension) == params.precision
    constants_polynomials.append(constants_mini_polynomial)
    constants_extensions.append(constants_mini_extension)
  assert len(constants_extensions) == deg
  for extension in constants_extensions:
    assert len(extension) == params.precision 
  print(
      'Converted round constants into a polynomial and low-degree extended it')
  return constants_extensions, constants_polynomials

def construct_computation_polynomial(comp, params):
  """Constructs polynomial for the given computation."""
  # Interpolate the computational trace into a polynomial P,
  # with each step along a successive power of G1
  #########################################################
  print("type(comp.computational_trace[0][0])")
  print(type(comp.computational_trace[0][0]))
  #assert 0 == 1
  #########################################################
  computational_trace_polynomial = fft(
      comp.computational_trace, params.modulus, params.G1,
      inv=True, dims=comp.dims)
  assert len(computational_trace_polynomial) == comp.steps
  p_evaluations = fft(computational_trace_polynomial,
      params.modulus, params.G2, dims=comp.dims)
  assert len(p_evaluations) == comp.steps*params.extension_factor
  print(
      'Converted computational steps into a polynomial and low-degree extended it'
  )
  return p_evaluations

def construct_constraint_polynomial(comp, params,
    p_evaluations):
  """Construct the constraint polynomial for the given tape.

  This function constructs a constraint polynomial for the
  given computational tape. For now, this function only works
  with MiMC.
  """
  #f = comp.field
  deg = len(comp.constants[0])
  extensions, _ = construct_constants_polynomials(comp, params)

  # Create the composed polynomial such that
  # C(P(x), P(g1*x), K(x)) = P(g1*x) - step_fn(P(x), K(x))
  # here K(x) contains the constants.
  p_next_step_evals = [p_evaluations[(i + params.extension_factor) % params.precision] for i in range(params.precision)]
  # extensions[d] selects constants for the degree d term
  # extensions[d][i] selects the degree d term for i-th step
  # extensions[d][i][0] unpacks the output of fft() which adds an extra list
  step_p_evals = [comp.step_fn(
    comp.field, p_evaluations[i], [extensions[d][i][0] for d in range(deg)]) for i in range(params.precision)]
  #c_of_p_evals = [[p_next[dim] - step_p[dim] % params.modulus for dim in range(comp.dims)] for (p_next, step_p) in zip(p_next_step_evals, step_p_evals)]
  c_of_p_evals = [[p_next[dim] - step_p[dim] for dim in range(comp.dims)] for (p_next, step_p) in zip(p_next_step_evals, step_p_evals)]
  print('Computed C(P, K) polynomial')
  return c_of_p_evals

def construct_remainder_polynomial(comp, params, c_of_p_evaluations):
  """Computes the remainder polynomial for the STARK.
  
  Compute D(x) = C(P(x), P(g1*x), K(x)) / Z(x)
  Z(x) = (x^steps - 1) / (x - x_atlast_step)
  TODO(rbharath): I think this is supposed to equal 
  Z(x) = (x - 1)(x-2)...(x-(steps_1)). How are these equal?
  """
  f = comp.field
  z_num_evaluations = [
      params.xs[(i * comp.steps) % params.precision] - 1 for i in range(params.precision)
  ]
  #z_num_inv = f.multi_inv(z_num_evaluations)
  z_num_inv = multi_inv(z_num_evaluations)
  # (x_i - x_{step-1}) list
  z_den_evaluations = [params.xs[i] - params.last_step_position for i in range(params.precision)]
  d_evaluations = [
      #[int(cp[dim] * zd * zni % params.modulus) for dim in range(comp.dims)]
      [cp[dim] * zd * zni for dim in range(comp.dims)]
      for cp, zd, zni in zip(c_of_p_evaluations, z_den_evaluations, z_num_inv)
  ]
  print('Computed D polynomial')
  return d_evaluations

def construct_boundary_polynomial(comp, params, p_evaluations):
  """Polynomial encoding boundary constraints on tape.
  
  Compute interpolant of ((1, input), (x_atlast_step, output))
  """
  field = comp.field
  polysOver = polynomials_over(field).factory
  i_evaluations = []
  inv_z2_evaluations = []
  #zeropoly2 = f.mul_polys([-1, 1], [-params.last_step_position, 1])
  zeropoly2 = polysOver([-1, 1])*polysOver([-params.last_step_position, 1])
  for dim in range(comp.dims):
    #interpolant = f.lagrange_interp_2([1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
    interpolant = lagrange_interp_2(params.modulus, polysOver([1, params.last_step_position]),
        polysOver([comp.inp[dim], comp.output[dim]]))
    #i_evaluations_dim = [f.eval_poly_at(interpolant, x) for x in params.xs]
    i_evaluations_dim = [interpolant(x) for x in params.xs]
    #inv_z2_evaluations_dim = f.multi_inv([f.eval_poly_at(zeropoly2, x) for x in params.xs])
    inv_z2_evaluations_dim = multi_inv([zeropoly2(x) for x in params.xs])
    # Append to list
    i_evaluations.append(i_evaluations_dim)
    inv_z2_evaluations.append(inv_z2_evaluations_dim)
  i_evaluations = [[i_evaluations[dim][j] for dim in range(comp.dims)] for j in range(params.precision)]
  inv_z2_evaluations = [[inv_z2_evaluations[dim][j] for dim in range(comp.dims)] for j in range(params.precision)]
  # B = (P - I) / Z2
  b_evaluations = []
  for p, i, invq in zip(p_evaluations, i_evaluations, inv_z2_evaluations):
    b_evaluations_dim = [
      #((p[dim] - i[dim]) * invq[dim]) % params.modulus for dim in range(comp.dims) ]
      (p[dim] - i[dim]) * invq[dim] for dim in range(comp.dims) ]
    b_evaluations.append(b_evaluations_dim)
  print('Computed B polynomial')
  return b_evaluations

def get_pseudorandom_ks(m_root, num):
  """Computes pseudorandom values from mtree root for linear combo.
  
  Parameters
  ----------
  m_root: bytestring
    Merkle tree root
  num: Int
    Number of pseudorandom positions to pull out
  """
  # Handling this case separately for backwards consistency. Might factor
  # out once more unit tests are in place.
  if 0 <= num and num <= 4:
    byte_list = [b'0x01', b'0x02', b'0x03', b'0x04']
    ks = [int.from_bytes(blake(m_root + byte_list[ind]), 'big') for ind in range(num)]
    return ks
  # Is the < 10 restriction needed
  elif num < 10:
    byte_list = [("0x0%s" % str(i)).encode("UTF-8") for i in range(num)]
    ks = [int.from_bytes(blake(m_root + byte_list[ind]), 'big') for ind in range(comp.dims)]
    return ks

def compute_pseudorandom_linear_combination_1d(comp, params, mtree, polys):
  """Computes the pseudorandom linear combination for 1-d slice of poly.

  A FRI proofs of low degree for a polynomial takes space.
  There are multiple polynomials (constraint, tape, reminder,
  degree) used in a STARK. This function combines these into a
  single polynomial so that only one FRI proofs needs to be
  generated for all of them. The chances of a collision are
  low.
  """
  f = comp.field
  # Based on the hashes of P, D and B, we select a random
  # linear combination of P * x^steps, P, B * x^steps, B and
  # D, and prove the low-degreeness of that, instead of
  # proving the low-degreeness of P, B and D separately
  k1, k2, k3, k4 = get_pseudorandom_ks(mtree[1], 4)
  # TODO(rbharath): This isn't general, but fix later
  [p_evaluations, d_evaluations, b_evaluations] = polys
  # Compute the linear combination. We don't even both
  # calculating it in coefficient form; we just compute the
  # evaluations
  #G2_to_the_steps = f.exp(params.G2, comp.steps)
  G2_to_the_steps = params.G2**comp.steps
  powers = [1]
  for i in range(1, params.precision):
    #powers.append(powers[-1] * G2_to_the_steps % params.modulus)
    powers.append(powers[-1] * G2_to_the_steps)

  l_evaluations_per_dim = []
  for dim in range(comp.dims):
    #l_evaluations_dim = [(d_evaluations[i][dim] + p_evaluations[i][dim] * k1 + p_evaluations[i][dim] * k2 * powers[i] + b_evaluations[i][dim] * k3 + b_evaluations[i][dim] * k4 * powers[i]) % params.modulus for i in range(params.precision)]
    l_evaluations_dim = [(d_evaluations[i][dim] + p_evaluations[i][dim] * k1 + p_evaluations[i][dim] * k2 * powers[i] + b_evaluations[i][dim] * k3 + b_evaluations[i][dim] * k4 * powers[i]) for i in range(params.precision)]
    l_evaluations_per_dim.append(l_evaluations_dim)
  return l_evaluations_per_dim

def compute_pseudorandom_linear_combination(comp, params, mtree, polys):
  """Computes a pseudorandom linear combination of polys

  A deterministic procedure for pseudorandomly combining dimensions
  """
  f = comp.field
  #G2_to_the_steps = f.exp(params.G2, comp.steps)
  G2_to_the_steps = params.G2**comp.steps
  powers = [1]
  for i in range(1, params.precision):
    #powers.append(powers[-1] * G2_to_the_steps % params.modulus)
    powers.append(powers[-1] * G2_to_the_steps)
  l_evaluations_per_dim = compute_pseudorandom_linear_combination_1d(comp,
      params, mtree, polys)
  l_ks = get_pseudorandom_ks(mtree[1], comp.dims)
  #l_evaluations = [sum([l_evals_dim[i] + l_evals_dim[i] * l_k * powers[i] for (l_evals_dim, l_k) in zip(l_evaluations_per_dim, l_ks)]) % params.modulus for i in range(params.precision)]
  l_evaluations = [sum([l_evals_dim[i] + l_evals_dim[i] * l_k * powers[i] for (l_evals_dim, l_k) in zip(l_evaluations_per_dim, l_ks)]) for i in range(params.precision)]
  print('Computed random linear combination')
  return l_evaluations

# TODO(rbharath): This function is poorly structured since it computes spot
# checks for both the mtree and the ltree simultaneously. This makes
# refactoring challenging. Break up and separate in future PR.
def compute_merkle_spot_checks(mtree, l_mtree, comp, params, samples=80):
  """Computes pseudorandom spot checks of Merkle tree."""
  # Do some spot checks of the Merkle tree at pseudo-random
  # coordinates, excluding multiples of `extension_factor`
  branches = []
  positions = get_pseudorandom_indices(
      l_mtree[1], params.precision, samples, exclude_multiples_of=params.extension_factor)
  for pos in positions:
    branches.append(mk_branch(mtree, pos))
    branches.append(mk_branch(mtree, (pos + params.extension_factor) % params.precision))
    branches.append(mk_branch(l_mtree, pos))
  print('Computed %d spot checks' % samples)
  return branches


def mk_proof(comp, params):
  """Generate a STARK for a MIMC calculation
  
  Parameters
  ----------
  comp: A Computation object
    Stores the Algebraic Intermediate Representation
  modulus: Int
    TODO(rbharath): This shouldn't be exposed
  """
  start_time = time.time()

  p_evaluations = construct_computation_polynomial(
      comp, params)

  # Construct the constraint polynomial (represented as a list
  # of point evaluations)
  c_of_p_evaluations = construct_constraint_polynomial(
      comp, params, p_evaluations)

  d_evaluations = construct_remainder_polynomial(
      comp, params, c_of_p_evaluations)

  b_evaluations = construct_boundary_polynomial(
      comp, params, p_evaluations)

  polys = [p_evaluations, d_evaluations, b_evaluations]
  # Compute their Merkle root
  mtree = merkelize_polynomials(comp.dims, polys)

  l_evaluations = compute_pseudorandom_linear_combination(
      comp, params, mtree, polys)
  l_mtree = merkelize(l_evaluations)

  branches = compute_merkle_spot_checks(mtree, l_mtree, comp, params)

  # Return the Merkle roots of P and D, the spot check Merkle
  # proofs, and low-degree proofs of P and D
  o = [
      mtree[1], l_mtree[1], branches,
      prove_low_degree(
          l_evaluations,
          params.G2,
          comp.steps * comp.constraint_degree,
          params.modulus,
          exclude_multiples_of=comp.extension_factor)
  ]
  print("STARK computed in %.4f sec" % (time.time() - start_time))
  return o

def verify_proof(comp, params, proof):
  """Verifies a STARK
  
  Parameters
  ----------
  comp: Computation 
    An Algebraic Intermediate Representation
    TODO(rbharath): This function should not see comp! This wouldn't be present
    in the real protocol.
  params: StarkParams
    TODO(rbharath): Is this needed?
  """
  start_time = time.time()
  m_root, l_root, branches, fri_proof = proof

  _, constants_polynomials = construct_constants_polynomials(comp, params)

  # Verifies the low-degree proofs
  assert verify_low_degree_proof(
      l_root,
      params.G2,
      fri_proof,
      comp.steps * comp.constraint_degree,
      params.modulus,
      exclude_multiples_of=comp.extension_factor)

  ## Performs the spot checks
  samples = params.spot_check_security_factor
  positions = get_pseudorandom_indices(
      l_root, params.precision, samples,
      exclude_multiples_of=params.extension_factor)
  ks = get_pseudorandom_ks(m_root, 4)
  for i, pos in enumerate(positions):
    verify_proof_at_position(comp, params, ks, proof, i, pos, constants_polynomials)

  print('Verified %d consistency checks' % params.spot_check_security_factor)
  print('Verified STARK in %.4f sec' % (time.time() - start_time))
  return True

def verify_proof_at_position(comp, params, ks, proof, i, pos, constants_polynomials):
  """Verifies merkle proof at given position in extended trace"""
  k1, k2, k3, k4 = ks
  m_root, l_root, branches, fri_proof = proof
  #x = comp.field.exp(params.G2, pos)
  x = params.G2**pos
  #x_to_the_steps = comp.field.exp(x, comp.steps)
  x_to_the_steps = x**comp.steps
  # Recall m is the merkle tree of the raw polynomials, and l
  # is the merkle tree of the pseudorandom combination
  # polynomial. Leaf node from m[pos]
  mbranch1 = verify_branch(m_root, pos, branches[i * 3])
  unpacked_leaf1 = unpack_merkle_leaf(mbranch1, comp.dims, 3)
  # Leaf node from m[pos + extension_factor]
  mbranch2 = verify_branch(
      m_root,
      (pos + params.extension_factor) % params.precision,
      branches[i * 3 + 1])
  unpacked_leaf2 = unpack_merkle_leaf(mbranch2, comp.dims, 3)
  # Leaf node from l[pos]
  l_of_x = verify_branch(l_root, pos, branches[i * 3 + 2],
      output_as_int=True)

  # This undoes the packing that's done in merkelize_polynomials
  # polys = [p_evaluations, d_evaluations, b_evaluations]
  p_of_x = [int.from_bytes(p_of_x_dim, 'big') for p_of_x_dim in unpacked_leaf1[:comp.dims]]
  p_of_g1x = [int.from_bytes(p_of_g1x_dim, 'big') for p_of_g1x_dim in unpacked_leaf2[:comp.dims]]
  d_of_x = [int.from_bytes(d_of_x_dim, 'big') for d_of_x_dim in unpacked_leaf1[comp.dims:2*comp.dims]]
  b_of_x = [int.from_bytes(b_of_x_dim, 'big') for b_of_x_dim in unpacked_leaf1[2*comp.dims:]]

  #zvalue = comp.field.div(comp.field.exp(x, comp.steps) - 1, x - params.last_step_position)
  zvalue = (x**comp.steps - 1)/(x - params.last_step_position)
  k_of_xs = []
  for constants_mini_polynomial in constants_polynomials:
    # This is unwrapping the polynomial
    constants_mini_polynomial = [val[0] for val in constants_mini_polynomial]
    k_of_x = comp.field.eval_poly_at(constants_mini_polynomial, x)
    k_of_xs.append(k_of_x)

  # Check transition constraints C(P(x)) = Z(x) * D(x)
  f_of_p_of_x = comp.step_fn(comp.field, p_of_x, k_of_xs)
  for dim in range(comp.dims):
    p_of_g1x_dim = p_of_g1x[dim]
    p_of_x_dim = p_of_x[dim]
    d_of_x_dim = d_of_x[dim]
    f_of_p_of_x_dim = f_of_p_of_x[dim]
    #assert (p_of_g1x_dim - f_of_p_of_x_dim - zvalue * d_of_x_dim) % params.modulus == 0
    assert (p_of_g1x_dim - f_of_p_of_x_dim - zvalue * d_of_x_dim) == 0

  # Check boundary constraints B(x) * Q(x) + I(x) = P(x)
  zeropoly2 = comp.field.mul_polys([-1, 1], [-params.last_step_position, 1])
  for dim in range(comp.dims):
    interpolant_dim = comp.field.lagrange_interp_2([1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
    #assert (p_of_x[dim] - b_of_x[dim] * comp.field.eval_poly_at(zeropoly2, x) - comp.field.eval_poly_at(interpolant_dim, x)) % params.modulus == 0
    assert (p_of_x[dim] - b_of_x[dim] * comp.field.eval_poly_at(zeropoly2, x) - comp.field.eval_poly_at(interpolant_dim, x)) == 0

  # TODO(rbharath): I'm commenting this out for now, but I think commenting
  # out this check breaks security guarantees!! To fix this, we need a way
  # of getting the dimensionwise l's, which might necessitate passing more
  # merkle branches into the original proof. Will refactor in a subsequent
  # PR.
  # Check correctness of the linear combination
  #assert (l_of_x - d_of_x - k1 * p_of_x - k2 * p_of_x * x_to_the_steps -
  #        k3 * b_of_x - k4 * b_of_x * x_to_the_steps) % modulus == 0
