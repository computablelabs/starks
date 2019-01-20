import time
from typing import List
from typing import Tuple
from starks.merkle_tree import blake
from starks.merkle_tree import verify_branch
from starks.merkle_tree import mk_branch
from starks.merkle_tree import merkelize
from starks.merkle_tree import merkelize_polynomials
from starks.merkle_tree import unpack_merkle_leaf
from starks.polynomial import polynomials_over
from starks.poly_utils import lagrange_interp_2 
from starks.fft import NonBinaryFFT
# TODO(rbharath): Swap these out with object oriented API 
#from starks.fri import prove_low_degree, verify_low_degree_proof
from starks.utils import generate_Xi_s
from starks.utils import get_power_cycle
from starks.utils import is_a_power_of_2
from starks.utils import get_pseudorandom_field_elements
from starks.air import Computation
from starks.poly_utils import make_multivar
from starks.poly_utils import multi_inv
from starks.numbertype import Field
from starks.numbertype import FieldElement
from starks.numbertype import Vector
from starks.numbertype import Poly

class StarkParams(object):
  """Holds the cryptographic parameters needed for STARK"""
  def __init__(self, field, steps: int, modulus: FieldElement,
      extension_factor: int, width: int, step_polys: List[Poly], spot_check_security_factor: int =80):
    """
    TODO(rbharath): I believe what this class is doing is
    constructing a smooth multiplicative group. Alternatively,
    this could be an affine space.

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
    self.width = width
    self.steps = steps
    self.step_polys = step_polys
    self.extension_factor = extension_factor
    self.precision = steps * extension_factor
    self.spot_check_security_factor = spot_check_security_factor

    # TODO(rbharath): Perhaps these should be roots of the primitive
    # polynomials in the full-fledged starks.
    # Root of unity such that x^precision=1
    self.G2 = field(7)**((modulus - 1) // self.precision)

    # Root of unity such that x^steps=1
    self.G1 = self.G2**extension_factor

    ## Powers of the higher-order root of unity
    #self.xs = get_power_cycle(self.G2, modulus)
    self.xs = get_power_cycle(self.G2, self.field)
    self.last_step_position = self.xs[(steps - 1) * extension_factor]

# TODO(rbharath): This is where the binary field fft will come into play.
def construct_trace_polynomials(witness, params: StarkParams) -> Poly:
  """Constructs polynomial for the given computation."""
  # Interpolate the computational trace into a polynomial P,
  # with each step along a successive power of G1
  # TODO(rbharath): This could be better...
  field = params.field
  nonbinary_fft = NonBinaryFFT(field, params.G1, params.width)
  trace_polys = []
  for witness_dim in witness:
    dim_trace = nonbinary_fft.inv_fft(witness_dim)
    trace_polys.append(dim_trace)
  return trace_polys

def construct_constraint_polynomials(trace_polys: List[Poly], params: StarkParams) -> List[Poly]:
  """Construct the constraint polynomial for the given tape.

  This function constructs a constraint polynomial for the
  given computational tape. For now, this function only works
  with MiMC.
  """
  # Create the composed polynomial such that
  # C(P(x), P(g1*x)) = P(g1*x) - step_fn(P(x))
  field, width = params.field, params.width
  Xi_s = generate_Xi_s(field, width)
  next_traces = [trace_poly(params.G1*X_i) for (trace_poly, X_i) in zip(trace_polys, Xi_s)]
  # Convert trace polys to multidimensional polys by evaluating
  trace_polys = [trace_poly(X_i) for (X_i, trace_poly) in zip(Xi_s, trace_polys)]
  constraint_polys = []
  transition_poly = step_poly(trace_polys)
  for next_trace, step_poly in zip(next_traces, params.step_polys):
    # TODO(rbharath): One of these is 1-d poly, while other is multi-d. How to handle?
    constraint_poly = next_trace - transition_poly
    constraint_polys.append(constraint_poly)
  return constraint_polys

def construct_remainder_polynomials(constraint_polys: List[Poly], params: StarkParams) -> List[Poly]:
  """Computes the remainder polynomial for the STARK.
  
  Compute D(x) = C(P(x), P(g1*x)) / Z(x)
  Z(x) = (x^steps - 1) / (x - x_atlast_step)
  TODO(rbharath): I think this is supposed to equal 
  Z(x) = (x - 1)(x-2)...(x-(steps_1)). How are these equal?
  """
  Xi_s = generate_Xi_s(params.field, params.width)
  polysOver = polynomials_over(params.field).factory
  X = polysOver([0, 1])
  # TODO(rbharath): Write a unit test checking that the behavior of z is as desired (x-1)...(x-(steps-1))
  z_num = X**params.steps - 1
  z_den = X - params.last_step_position
  # Check division is OK
  assert z_num % z_den == 0
  z = z_num / z_den
  # Implicitly representing the division...
  ds = [(cp, z(X_i)) for (cp, X_i) in zip(constraint_polys, Xi_s)]
  print('Computed D polynomials')
  #return d_evaluations
  return ds

#def construct_boundary_polynomials(comp: Computation, params: StarkParams, p_evaluations: List[Vector]) -> List[Vector]:
def construct_boundary_polynomials(trace_polys: List[Poly], witness: List[List], boundary: List[Tuple], params: StarkParams) -> List[Vector]:
  """Polynomial encoding boundary constraints on tape.
  
  Compute interpolant of ((1, input), (x_atlast_step, output))

  TODO(rbharath): This assumes boundary has simplified form
  """
  polysOver = polynomials_over(params.field).factory
  interpolants = []
  inv_z2_polys = []
  zeropoly2 = polysOver([-1, 1])*polysOver([-params.last_step_position, 1])
  for dim in range(params.width):
    constraint = boundary[dim]
    (_, _, input_value) = constraint
    output_dim = witness[dim][-1]
    interpolant = lagrange_interp_2(params.field, polysOver([1, params.last_step_position]),
        polysOver([input_value, output_dim]))
    interpolants.append(interpolant)
  # B = (P - I) / Z2
  b_polys = []
  for p, i in zip(trace_polys, interpolants):
    b_poly = [
      (p - i)/zeropoly2 for dim in range(params.width) ]
    #b_evaluations.append(b_evaluations_dim)
    b_polys.append(b_poly)
  print('Computed B polynomial')
  #return b_evaluations
  return b_polys

def get_pseudorandom_ks(m_root: bytes, num: int) -> List[int]:
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
    ks = [int.from_bytes(blake(m_root + byte_list[ind]), 'big') for ind in range(num)]
    return ks

def compute_pseudorandom_linear_combination_1d(params: StarkParams, mtree: List[bytes], polys: List[List[Vector]]) -> List[Vector]:
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
  #[p_evaluations, d_evaluations, b_evaluations] = polys
  [trace_polys, reminder_polys, boundary_polys] = polys
  # Compute the linear combination. We don't even both
  # calculating it in coefficient form; we just compute the
  # evaluations
  G2_to_the_steps = params.G2**params.steps
  powers = [1]
  for i in range(1, params.precision):
    powers.append(powers[-1] * G2_to_the_steps)

  #l_evaluations_per_dim = []
  l_polys = []
  for dim in range(params.width):
    #l_evaluations_dim = [(d_evaluations[i][dim] + p_evaluations[i][dim] * k1 + p_evaluations[i][dim] * k2 * powers[i] + b_evaluations[i][dim] * k3 + b_evaluations[i][dim] * k4 * powers[i]) for i in range(params.precision)]
    l_poly = reminder_poly + trace_poly * k1 + trace_poly * k2 * powers[i] + boundary_poly * k3 + boundary_poly * k4 * powers[i]
    #l_evaluations_per_dim.append(l_evaluations_dim)
    l_polys.append(l_poly)
  #return l_evaluations_per_dim
  return l_polys

def compute_pseudorandom_linear_combination(params: StarkParams, mtree: List[bytes], polys: List[List[Vector]]):
  """Computes a pseudorandom linear combination of polys

  A deterministic procedure for pseudorandomly combining dimensions
  """
  f = params.field
  G2_to_the_steps = params.G2**params.steps
  powers = [1]
  for i in range(1, params.precision):
    powers.append(powers[-1] * G2_to_the_steps)
  l_evaluations_per_dim = compute_pseudorandom_linear_combination_1d(params, mtree, polys)
  l_ks = get_pseudorandom_ks(mtree[1], params.width)
  l_evaluations = [sum([l_evals_dim[i] + l_evals_dim[i] * l_k * powers[i] for (l_evals_dim, l_k) in zip(l_evaluations_per_dim, l_ks)]) for i in range(params.precision)]
  print('Computed random linear combination')
  return l_evaluations

# TODO(rbharath): This function is poorly structured since it computes spot
# checks for both the mtree and the ltree simultaneously. This makes
# refactoring challenging. Break up and separate in future PR.
def compute_merkle_spot_checks(mtree, l_mtree, params, samples=80):
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


#def mk_proof(comp: Computation, params: StarkParams):
def mk_proof(comp: Computation, params: StarkParams):
  """Generate a STARK for a MIMC calculation"""
  start_time = time.time()

  fft_solver = FFT(params.field)
  # TODO(rbharath): I think by computing with explicit
  # polynomials, this function becomes a *lot* simpler. Then
  # the evaluation can be done only once.
  #p_evaluations = construct_computation_polynomial(
  #    comp, params)

  # Construct the constraint polynomial (represented as a list
  # of point evaluations)
  #c_of_p_evaluations = construct_constraint_polynomial(
  #    comp, params, p_evaluations)

  # list of length |width|
  trace_polys = construct_trace_polys(witness, params)
  constraint_polys = construct_constraint_polynomials(
      trace_polys, params)
  remainder_polys = construct_remainder_polynomials(
      constraint_polys, params)
  boundary_polys = construct_boundary_polynomials(
      comp, params, p_evaluations)

  #polys = [p_evaluations, d_evaluations, b_evaluations]
  polys = [trace_polys, remainder_polys, boundary_polys]
  # Compute their Merkle root
  # TODO(rbharath): The merkelization is computed on the
  # affine subspace of the RS[F, L, pho] I believe.
  # Alternatively on a smooth multiplicative group, which is
  # what's happening now.
  mtree = merkelize_polynomials(comp.width, polys)

  #l_evaluations = compute_pseudorandom_linear_combination(
  #    comp, params, mtree, polys)
  l_polys = compute_pseudorandom_linear_combination(
      comp, params, mtree, polys)
  l_mtree = merkelize(l_evaluations)

  branches = compute_merkle_spot_checks(mtree, l_mtree, comp, params)

  fri = FRI()
  # Return the Merkle roots of P and D, the spot check Merkle
  # proofs, and low-degree proofs of P and D
  o = [
      mtree[1], l_mtree[1], branches,
      #prove_low_degree(
      #    #l_evaluations,
      #    l_polys,
      #    params.G2,
      #    comp.steps * comp.get_degree(),
      #    #params.modulus,
      #    params.field, # TODO(rbharath): This should be serialized
      #    exclude_multiples_of=comp.extension_factor)
      fri.prove_proximity(l_polys, params.G2, comp.steps*comp.get_degree())
  ]
  print("STARK computed in %.4f sec" % (time.time() - start_time))
  return o

def verify_proof(comp: Computation, params: StarkParams, proof):
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

  #_, constants_polynomials = construct_constants_polynomials(comp, params)

  # Verifies the low-degree proofs
  assert verify_low_degree_proof(
      l_root,
      params.G2,
      fri_proof,
      comp.steps * comp.get_degree(),
      #params.modulus,
      params.field,
      exclude_multiples_of=comp.extension_factor)

  ## Performs the spot checks
  samples = params.spot_check_security_factor
  positions = get_pseudorandom_indices(
      l_root, params.precision, samples,
      exclude_multiples_of=params.extension_factor)
  ks = get_pseudorandom_ks(m_root, 4)
  for i, pos in enumerate(positions):
    #verify_proof_at_position(comp, params, ks, proof, i, pos, constants_polynomials)
    verify_proof_at_position(comp, params, ks, proof, i, pos)

  print('Verified %d consistency checks' % params.spot_check_security_factor)
  print('Verified STARK in %.4f sec' % (time.time() - start_time))
  return True

def verify_proof_at_position(comp, params, ks, proof, i, pos):
  """Verifies merkle proof at given position in extended trace"""
  field = comp.field
  #modulus = params.modulus
  polysOver = polynomials_over(field).factory
  k1, k2, k3, k4 = ks
  m_root, l_root, branches, fri_proof = proof
  x = params.G2**pos
  x_to_the_steps = x**comp.steps
  # Recall m is the merkle tree of the raw polynomials, and l
  # is the merkle tree of the pseudorandom combination
  # polynomial. Leaf node from m[pos]
  mbranch1 = verify_branch(m_root, pos, branches[i * 3])
  unpacked_leaf1 = unpack_merkle_leaf(mbranch1, comp.width, 3)
  # Leaf node from m[pos + extension_factor]
  mbranch2 = verify_branch(
      m_root,
      (pos + params.extension_factor) % params.precision,
      branches[i * 3 + 1])
  unpacked_leaf2 = unpack_merkle_leaf(mbranch2, comp.width, 3)
  # Leaf node from l[pos]
  l_of_x = verify_branch(l_root, pos, branches[i * 3 + 2],
      output_as_int=True)

  # This undoes the packing that's done in merkelize_polynomials
  p_of_x = [field(p_of_x_dim) for p_of_x_dim in unpacked_leaf1[:comp.width]]
  p_of_g1x = [field(p_of_g1x_dim) for p_of_g1x_dim in unpacked_leaf2[:comp.width]]
  d_of_x = [field(d_of_x_dim) for d_of_x_dim in unpacked_leaf1[comp.width:2*comp.width]]
  b_of_x = [field(b_of_x_dim) for b_of_x_dim in unpacked_leaf1[2*comp.width:]]

  zvalue = (x**comp.steps - 1)/(x - params.last_step_position)
  k_of_xs = []

  # Check transition constraints C(P(x)) = Z(x) * D(x)
  f_of_p_of_x = [comp.step_polys[i](p_of_x) for i in range(comp.width)]
  for dim in range(comp.width):
    p_of_g1x_dim = p_of_g1x[dim]
    p_of_x_dim = p_of_x[dim]
    d_of_x_dim = d_of_x[dim]
    f_of_p_of_x_dim = f_of_p_of_x[dim]
    assert (p_of_g1x_dim - f_of_p_of_x_dim - zvalue * d_of_x_dim) == 0

  # Check boundary constraints B(x) * Q(x) + I(x) = P(x)
  # TODO(rbharath): How do I promote a single-dim poly into a multidimensional poly?
  zeropoly2 = polysOver([-1, 1])*polysOver([-params.last_step_position, 1])
  for dim in range(comp.width):
    #interpolant_dim = lagrange_interp_2(modulus, [1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
    interpolant_dim = lagrange_interp_2(field, [1, params.last_step_position], [comp.inp[dim], comp.output[dim]])
    assert (p_of_x[dim] - b_of_x[dim] * zeropoly2(x) - interpolant_dim(x)) == 0

  # TODO(rbharath): I'm commenting this out for now, but I think commenting
  # out this check breaks security guarantees!! To fix this, we need a way
  # of getting the dimensionwise l's, which might necessitate passing more
  # merkle branches into the original proof. Will refactor in a subsequent
  # PR.
  # Check correctness of the linear combination
  #assert (l_of_x - d_of_x - k1 * p_of_x - k2 * p_of_x * x_to_the_steps -
  #        k3 * b_of_x - k4 * b_of_x * x_to_the_steps) % modulus == 0
