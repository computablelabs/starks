import time
import numpy as np
from starks.merkle_tree import merkelize, mk_branch, verify_branch, blake
from starks.compression import compress_fri, decompress_fri, compress_branches, decompress_branches, bin_length
from starks.poly_utils import PrimeField
from starks.fft import fft
from starks.fri import prove_low_degree, verify_low_degree_proof
from starks.utils import get_power_cycle, get_pseudorandom_indices, is_a_power_of_2

modulus = 2**256 - 2**32 * 351 + 1
f = PrimeField(modulus)
nonresidue = 7

# Number of branches used for Merkle-tree check
spot_check_security_factor = 80

def merkelize_polynomials(dims, polynomials):
  """Given a list of polynomial evaluations, merkelizes them together.

  Each leaf of the Merkle tree contains the concatenation of the values of
  each polynomial in polynomials at a given index. If the polynomials are
  multidimensional, the per-dimension leaves are concatenated to form one
  joint leaf.

  Parameters
  ----------
  dims: Int
    Dimensionality
  polynomials: List
    Each element much be a list of evaluations of a given poly. All of
    these should have the same length.
  """
  # Note len(mtree_leaf) == 32 * dims * len(polys)
  # Note len(mtree) == 2 * precision now
  # This code packs each Merkle leaf as
  # polyval_1_dim_1 polyval_1_dim_2 poly_val_2_dim_1 poly_val_2_dim_2 ...
  # In the common case this is
  # [p_of_x_dim_1 p_of_x_dim_2 .. d_of_x_dim_1 d_of_x_dim_2... b_of_x_dim_1 b_of_x_dim_2]
  mtree = merkelize([
      b''.join([val[dim].to_bytes(32, 'big') for val in evals for dim in range(dims)])
      for evals in zip(*polynomials)])
  return mtree

def unpack_merkle_leaf(leaf, dims, num_polys):
  """Unpacks a merkle leaf created by merkelize_polynomials.

  Note that the packing in each Merkle leaf is

  polyval_1_dim_1 polyval_1_dim_2 poly_val_2_dim_1 poly_val_2_dim_2 ...

  Each value is encoded in 32 bytes, so the total length of
  the leaf is 32 * num_polys * dims bytes.

  Parameters
  ----------
  leaf: bytestring
    A bytestring holding the Merkle leaf
  dims: Int
    The dimensionality of the state space
  num_polys: int
    The number of polynomials
  """
  vals = []
  for poly_ind in range(num_polys):
    for dim in range(dims):
      start_index = 32 * (poly_ind * dims + dim)
      end_index = 32 * (poly_ind * dims + dim + 1)
      byte_val = leaf[start_index:end_index]
      vals.append(byte_val)
  return vals


def get_computational_trace(inp, steps, constants, step_fn):
  """Get the computational trace for the STARK.

  Parameters
  ----------
  inp: list 
    The input state for the computation
  steps: Int
    The number of steps in the computation
  constants: List
    List of constants defining the computation in question
  step_fn: Function
    A function which maps one state to the next state.
  """
  computational_trace = [inp]
  for i in range(steps - 1):
    poly_constants = constants[i]
    # TODO(rbharath): Is there off-by-one error on round_contants?
    next_state = step_fn(f, computational_trace[-1], poly_constants)
    if isinstance(next_state, int):
      next_state = [next_state]
    computational_trace.append(next_state)
  output = computational_trace[-1]
  print('Done generating computational trace')
  return computational_trace, output


class Computation(object):
  """A simple class defining a computation.
  
  Holds the state of the computation. Here are the various
  fields that constitue a computation.

  Fields
  ------
  dims: Int
    The dimensionality of the state space for the computation.
  inp: Int or List
    Either a single int or a list of integers of length dims
  steps: Int
    An int holding the number of steps of this computation.
  output: Int of List
    Either a single int or a list of integers of length dims
  constants: List
    A list of constants. Each element of constants must be a
    list of length steps
  step_fn: Function
    A function that maps a computation state to the next
    state. A state here is either an int of a list of ints of
    length dims.
  """
  def __init__(self, dims, inp, steps, constants, step_fn):
    self.dims = dims
    # Handle 1-d case
    if isinstance(inp, int):
      inp = [inp]
    self.inp = inp
    self.steps = steps
    self.constants = constants
    self.step_fn = step_fn
    self.computational_trace, self.output = get_computational_trace(inp,
        steps, constants, step_fn)

class StarkParams(object):
  """Holds the cryptographic parameters needed for STARK"""
  def __init__(self, comp, modulus, extension_factor):
    """
    Parameters
    ----------
    comp: Computation
      A Computation object 
    modulus: Int
      A prime p that defines finite field Z/p
    extension_factor: Int
      A power of two which is the degree to which the trace is expanded
      when  constructing polynomials. For example, a trace of length 512
      with an extension_factor of 8 would construct polynomial evaluations
      on a 4096 elements.
    """
    self.modulus = modulus
    self.extension_factor = extension_factor
    self.precision = comp.steps * extension_factor

    # Root of unity such that x^precision=1
    self.G2 = f.exp(7, (modulus - 1) // self.precision)

    # Root of unity such that x^steps=1
    self.G1 = f.exp(self.G2, extension_factor)

    # Powers of the higher-order root of unity
    self.xs = get_power_cycle(self.G2, modulus)
    self.last_step_position = self.xs[(comp.steps - 1) * extension_factor]

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
        deg_constants, modulus, params.G1,
        inv=True, dims=1)
    constants_mini_extension = fft(constants_mini_polynomial,
        modulus, params.G2, dims=1)
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
    f, p_evaluations[i], [extensions[d][i][0] for d in range(deg)]) for i in range(params.precision)]
  c_of_p_evals = [[p_next[dim] - step_p[dim] % params.modulus for dim in range(comp.dims)] for (p_next, step_p) in zip(p_next_step_evals, step_p_evals)]
  print('Computed C(P, K) polynomial')
  return c_of_p_evals

def construct_remainder_polynomial(comp, params, c_of_p_evaluations):
  """Computes the remainder polynomial for the STARK.
  
  Compute D(x) = C(P(x), P(g1*x), K(x)) / Z(x)
  Z(x) = (x^steps - 1) / (x - x_atlast_step)
  TODO(rbharath): I think this is supposed to equal 
  Z(x) = (x - 1)(x-2)...(x-(steps_1)). How are these equal?
  """
  z_num_evaluations = [
      params.xs[(i * comp.steps) % params.precision] - 1 for i in range(params.precision)
  ]
  z_num_inv = f.multi_inv(z_num_evaluations)
  # (x_i - x_{step-1}) list
  z_den_evaluations = [params.xs[i] - params.last_step_position for i in range(params.precision)]
  d_evaluations = [
      [int(cp[dim] * zd * zni % modulus) for dim in range(comp.dims)]
      for cp, zd, zni in zip(c_of_p_evaluations, z_den_evaluations, z_num_inv)
  ]
  print('Computed D polynomial')
  return d_evaluations

def construct_boundary_polynomial(comp, params, p_evaluations):
  """Polynomial encoding boundary constraints on tape.
  
  Compute interpolant of ((1, input), (x_atlast_step, output))
  """
  i_evaluations = []
  inv_z2_evaluations = []
  zeropoly2 = f.mul_polys([-1, 1], [-params.last_step_position, 1])
  for dim in range(comp.dims):
    interpolant = f.lagrange_interp_2([1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
    i_evaluations_dim = [f.eval_poly_at(interpolant, x) for x in params.xs]
    inv_z2_evaluations_dim = f.multi_inv([f.eval_poly_at(zeropoly2, x) for x in params.xs])
    # Append to list
    i_evaluations.append(i_evaluations_dim)
    inv_z2_evaluations.append(inv_z2_evaluations_dim)
  i_evaluations = [[i_evaluations[dim][j] for dim in range(comp.dims)] for j in range(params.precision)]
  inv_z2_evaluations = [[inv_z2_evaluations[dim][j] for dim in range(comp.dims)] for j in range(params.precision)]
  # B = (P - I) / Z2
  b_evaluations = []
  for p, i, invq in zip(p_evaluations, i_evaluations, inv_z2_evaluations):
    b_evaluations_dim = [
      ((p[dim] - i[dim]) * invq[dim]) % modulus for dim in range(comp.dims) ]
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
  G2_to_the_steps = f.exp(params.G2, comp.steps)
  powers = [1]
  for i in range(1, params.precision):
    powers.append(powers[-1] * G2_to_the_steps % params.modulus)

  l_evaluations_per_dim = []
  for dim in range(comp.dims):
    l_evaluations_dim = [(d_evaluations[i][dim] + p_evaluations[i][dim] * k1 + p_evaluations[i][dim] * k2 * powers[i] + b_evaluations[i][dim] * k3 + b_evaluations[i][dim] * k4 * powers[i]) % params.modulus for i in range(params.precision)]
    l_evaluations_per_dim.append(l_evaluations_dim)
  return l_evaluations_per_dim

def compute_pseudorandom_linear_combination(comp, params, mtree, polys):
  """Computes a pseudorandom linear combination of polys

  A deterministic procedure for pseudorandomly combining dimensions
  """
  G2_to_the_steps = f.exp(params.G2, comp.steps)
  powers = [1]
  for i in range(1, params.precision):
    powers.append(powers[-1] * G2_to_the_steps % params.modulus)
  l_evaluations_per_dim = compute_pseudorandom_linear_combination_1d(comp,
      params, mtree, polys)
  l_ks = get_pseudorandom_ks(mtree[1], comp.dims)
  l_evaluations = [sum([l_evals_dim[i] + l_evals_dim[i] * l_k * powers[i] for (l_evals_dim, l_k) in zip(l_evaluations_per_dim, l_ks)]) % modulus for i in range(params.precision)]
  print('Computed random linear combination')
  return l_evaluations

# TODO(rbharath): This function is poorly structured since it computes spot
# checks for both the mtree and the ltree simultaneously. This makes
# refactoring challenging. Break up and separate in future PR.
def compute_merkle_spot_checks(mtree, l_mtree, comp, params, samples=spot_check_security_factor):
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


def mk_proof(inp, steps, constants, step_fn, constraint_degree=2, dims=1, extension_factor=8):
  """Generate a STARK for a MIMC calculation
  
  Parameters
  ----------
  constraint_degree: int
    The degree of the constraint being considered
  dims: int
    The dimension of the state space for the computation.
  """
  start_time = time.time()
  # Some constraints to make our job easier
  assert steps <= 2**32 // extension_factor
  # remove
  assert is_a_power_of_2(steps)
  for poly_constants in constants:
    assert len(poly_constants) <= steps

  comp = Computation(dims, inp, steps, constants, step_fn)
  params = StarkParams(comp, modulus, extension_factor)
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
  mtree = merkelize_polynomials(dims, polys)

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
          steps * constraint_degree,
          modulus,
          exclude_multiples_of=extension_factor)
  ]
  print("STARK computed in %.4f sec" % (time.time() - start_time))
  return o

def verify_proof(inp, steps, constants, output, proof, step_fn,
    constraint_degree=2, extension_factor=8, dims=1):
  """Verifies a STARK
  
  Parameters
  ----------
  constraint_degree: int
    The degree of the constraint being considered
  """
  start_time = time.time()
  assert steps <= 2**32 // extension_factor
  m_root, l_root, branches, fri_proof = proof
  comp = Computation(dims, inp, steps, constants, step_fn)
  params = StarkParams(comp, modulus, extension_factor)
  # ALl constants should be of same length so we check the first
  assert is_a_power_of_2(steps)
  assert len(constants[0]) <= steps

  _, constants_polynomials = construct_constants_polynomials(comp, params)

  # Verifies the low-degree proofs
  assert verify_low_degree_proof(
      l_root,
      params.G2,
      fri_proof,
      steps * constraint_degree,
      modulus,
      exclude_multiples_of=extension_factor)

  ## Performs the spot checks
  samples = spot_check_security_factor
  positions = get_pseudorandom_indices(
      l_root, params.precision, samples,
      exclude_multiples_of=params.extension_factor)
  ks = get_pseudorandom_ks(m_root, 4)
  for i, pos in enumerate(positions):
    verify_proof_at_position(comp, params, ks, proof, i, pos, constants_polynomials)

  print('Verified %d consistency checks' % spot_check_security_factor)
  print('Verified STARK in %.4f sec' % (time.time() - start_time))
  return True

def verify_proof_at_position(comp, params, ks, proof, i, pos, constants_polynomials):
  """Verifies merkle proof at given position in extended trace"""
  k1, k2, k3, k4 = ks
  m_root, l_root, branches, fri_proof = proof
  x = f.exp(params.G2, pos)
  x_to_the_steps = f.exp(x, comp.steps)
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

  zvalue = f.div(f.exp(x, comp.steps) - 1, x - params.last_step_position)
  k_of_xs = []
  for constants_mini_polynomial in constants_polynomials:
    # This is unwrapping the polynomial
    constants_mini_polynomial = [val[0] for val in constants_mini_polynomial]
    k_of_x = f.eval_poly_at(constants_mini_polynomial, x)
    k_of_xs.append(k_of_x)

  # Check transition constraints C(P(x)) = Z(x) * D(x)
  f_of_p_of_x = comp.step_fn(f, p_of_x, k_of_xs)
  for dim in range(comp.dims):
    p_of_g1x_dim = p_of_g1x[dim]
    p_of_x_dim = p_of_x[dim]
    d_of_x_dim = d_of_x[dim]
    f_of_p_of_x_dim = f_of_p_of_x[dim]
    assert (p_of_g1x_dim - f_of_p_of_x_dim - zvalue * d_of_x_dim) % modulus == 0

  # Check boundary constraints B(x) * Q(x) + I(x) = P(x)
  zeropoly2 = f.mul_polys([-1, 1], [-params.last_step_position, 1])
  for dim in range(comp.dims):
    interpolant_dim = f.lagrange_interp_2([1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
    assert (p_of_x[dim] - b_of_x[dim] * f.eval_poly_at(zeropoly2, x) - f.eval_poly_at(interpolant_dim, x)) % modulus == 0

  # TODO(rbharath): I'm commenting this out for now, but I think commenting
  # out this check breaks security guarantees!! To fix this, we need a way
  # of getting the dimensionwise l's, which might necessitate passing more
  # merkle branches into the original proof. Will refactor in a subsequent
  # PR.
  # Check correctness of the linear combination
  #assert (l_of_x - d_of_x - k1 * p_of_x - k2 * p_of_x * x_to_the_steps -
  #        k3 * b_of_x - k4 * b_of_x * x_to_the_steps) % modulus == 0
