import time
from starks.merkle_tree import merkelize, mk_branch, verify_branch, blake
from starks.compression import compress_fri, decompress_fri, compress_branches, decompress_branches, bin_length
from starks.poly_utils import PrimeField
from starks.fft import fft
from starks.fri import prove_low_degree, verify_low_degree_proof
from starks.utils import get_power_cycle, get_pseudorandom_indices, is_a_power_of_2

modulus = 2**256 - 2**32 * 351 + 1
f = PrimeField(modulus)
nonresidue = 7

# Number of bits of security in Merkle-tree check
spot_check_security_factor = 80
# TODO(rbharath): Is this the Galois extension degree?
# I think this is a security factor that gives us more control
# over the security level of the computation, but not sure.
extension_factor = 8

def get_computational_trace(inp, steps, constants, computational_step):
  """Get the computational trace for the STARK.

  This function is a first step towards refactoring this code
  so it can generate STARKs for different computations. For
  now, this only works for MiMC

  Parameters
  ----------
  inp: Int
    The input for the computation
  """
  computational_trace = [inp]
  deg = len(constants)
  for i in range(steps - 1):
    poly_constants = [constants[d][i] for d in range(deg)]
    # TODO(rbharath): Is there off-by-one error on round_contants?
    computational_trace.append(computational_step(f, computational_trace[-1], poly_constants))
  output = computational_trace[-1]
  print('Done generating computational trace')
  return computational_trace, output

def construct_constraint_polynomial(steps, constants, G1, G2, precision, p_evaluations, step_fn):
  """Construct the constraint polynomial for the given tape.

  This function constructs a constraint polynomial for the
  given computational tape. For now, this function only works
  with MiMC.
  """
  deg = len(constants)
  constants_extensions = []
  for d in range(deg):
    deg_constants = constants[d]
    #deg_constants = constants[0]
    skips2 = steps // len(deg_constants)
    constants_mini_polynomial = fft(
        deg_constants, modulus, f.exp(G1, skips2), inv=True)
    constants_mini_extension = fft(constants_mini_polynomial, modulus,
                                  f.exp(G2, skips2))
    assert len(constants_mini_extension) == precision // skips2
    constants_extensions.append(constants_mini_extension)
  assert len(constants_extensions) == deg
  for extension in constants_extensions:
    assert len(extension) == precision // skips2
  print(
      'Converted round constants into a polynomial and low-degree extended it')

  # Create the composed polynomial such that
  # C(P(x), P(g1*x), K(x)) = P(g1*x) - P(x)**3 - K(x)
  # here K(x) is the round constants.
  c_of_p_evaluations = [
      (p_evaluations[(i + extension_factor) % precision]
        - step_fn(f, p_evaluations[i],
                  [constants_extensions[d][i] for d in range(deg)])
       ) % modulus
      for i in range(precision)
  ]
  print('Computed C(P, K) polynomial')
  return c_of_p_evaluations

def compute_remainder_polynomial(xs, precision, steps, last_step_position, c_of_p_evaluations):
  """Computes the remainder polynomial for the STARK.
  
  Compute D(x) = C(P(x), P(g1*x), K(x)) / Z(x)
  Z(x) = (x^steps - 1) / (x - x_atlast_step)
  TODO(rbharath): I think this is supposed to equal 
  Z(x) = (x - 1)(x-2)...(x-(steps_1)). How are these equal?
  """
  z_num_evaluations = [
      xs[(i * steps) % precision] - 1 for i in range(precision)
  ]
  z_num_inv = f.multi_inv(z_num_evaluations)
  # (x_i - x_{step-1}) list
  z_den_evaluations = [xs[i] - last_step_position for i in range(precision)]
  d_evaluations = [
      cp * zd * zni % modulus
      for cp, zd, zni in zip(c_of_p_evaluations, z_den_evaluations, z_num_inv)
  ]
  print('Computed D polynomial')
  return d_evaluations

def compute_boundary_polynomial(xs, last_step_position, inp, output, p_evaluations):
  """Polynomial encoding boundary constraints on tape.
  
  Compute interpolant of ((1, input), (x_atlast_step, output))
  """
  # TODO(rbharath): Why does this interpolant make sense?
  interpolant = f.lagrange_interp_2([1, last_step_position], [inp, output])
  i_evaluations = [f.eval_poly_at(interpolant, x) for x in xs]

  zeropoly2 = f.mul_polys([-1, 1], [-last_step_position, 1])
  inv_z2_evaluations = f.multi_inv([f.eval_poly_at(zeropoly2, x) for x in xs])

  # B = (P - I) / Z2
  b_evaluations = [
      ((p - i) * invq) % modulus
      for p, i, invq in zip(p_evaluations, i_evaluations, inv_z2_evaluations)
  ]
  print('Computed B polynomial')
  return b_evaluations


def compute_pseudorandom_linear_combination(mtree, G2, steps, precision, d_evaluations, p_evaluations, b_evaluations):
  """Computes a pseudorandom linear combination of polys

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
  k1 = int.from_bytes(blake(mtree[1] + b'\x01'), 'big')
  k2 = int.from_bytes(blake(mtree[1] + b'\x02'), 'big')
  k3 = int.from_bytes(blake(mtree[1] + b'\x03'), 'big')
  k4 = int.from_bytes(blake(mtree[1] + b'\x04'), 'big')

  # Compute the linear combination. We don't even both
  # calculating it in coefficient form; we just compute the
  # evaluations
  G2_to_the_steps = f.exp(G2, steps)
  powers = [1]
  for i in range(1, precision):
    powers.append(powers[-1] * G2_to_the_steps % modulus)

  l_evaluations = [(d_evaluations[i] + p_evaluations[i] * k1 +
                    p_evaluations[i] * k2 * powers[i] + b_evaluations[i] * k3 +
                    b_evaluations[i] * powers[i] * k4) % modulus
                   for i in range(precision)]

  print('Computed random linear combination')
  return l_evaluations

def compute_merkle_spot_checks(mtree, l_mtree, precision, skips):
  """Computes pseudorandom spot checks of Merkle tree."""
  # Do some spot checks of the Merkle tree at pseudo-random
  # coordinates, excluding multiples of `extension_factor`
  branches = []
  samples = spot_check_security_factor
  positions = get_pseudorandom_indices(
      l_mtree[1], precision, samples, exclude_multiples_of=extension_factor)
  for pos in positions:
    branches.append(mk_branch(mtree, pos))
    branches.append(mk_branch(mtree, (pos + skips) % precision))
    branches.append(mk_branch(l_mtree, pos))
  print('Computed %d spot checks' % samples)
  return branches


# NOTE(rbharath): If you run a deep learning model on a GPU,
# you can add a trace-log which can be exited from the GPU
# without much effort. This trace log can be passed along to
# the STARK prover off-line.

# TODO(rbharath): I think the general strategy is to create a
# "comptutational "trace" of the workload, then constract a
# polynomial constraint with necessary linkage. Easier for
# regular workloads. This is also called a "computation tape"

def mk_proof(inp, steps, constants, step_fn, constraint_degree=2, dims=1):
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
  # TODO(rbharath): The second check here may be feasible to
  # remove
  #assert is_a_power_of_2(steps) and is_a_power_of_2(len(constants))
  assert is_a_power_of_2(steps)
  for poly_constants in constants:
    assert len(poly_constants) <= steps

  precision = steps * extension_factor

  # Root of unity such that x^precision=1
  G2 = f.exp(7, (modulus - 1) // precision)

  # Root of unity such that x^steps=1
  skips = precision // steps
  assert skips == extension_factor
  G1 = f.exp(G2, skips)

  # Powers of the higher-order root of unity
  xs = get_power_cycle(G2, modulus)
  last_step_position = xs[(steps - 1) * extension_factor]

  # computational_trace is a tape of computation values. (Put
  # another way, a list of partial values the computation
  # takes on).
  computational_trace, output = get_computational_trace(inp, steps, constants, step_fn)

  # Interpolate the computational trace into a polynomial P,
  # with each step along a successive power of G1
  computational_trace_polynomial = fft(
      computational_trace, modulus, G1, inv=True, dims=dims)
  assert len(computational_trace_polynomial) == steps
  p_evaluations = fft(computational_trace_polynomial, modulus, G2, dims=dims)
  assert len(p_evaluations) == steps*extension_factor
  print(
      'Converted computational steps into a polynomial and low-degree extended it'
  )

  # Construct the constraint polynomial (represented as a list
  # of point evaluations)
  c_of_p_evaluations = construct_constraint_polynomial(steps,
      constants, G1, G2, precision, p_evaluations,
      step_fn)

  d_evaluations = compute_remainder_polynomial(xs, precision,
      steps, last_step_position, c_of_p_evaluations)

  b_evaluations = compute_boundary_polynomial(xs,
      last_step_position, inp, output, p_evaluations)

  # Compute their Merkle root
  mtree = merkelize([
      pval.to_bytes(32, 'big') + dval.to_bytes(32, 'big') + bval.to_bytes(
          32, 'big')
      for pval, dval, bval in zip(p_evaluations, d_evaluations, b_evaluations)
  ])
  print('Computed hash root')

  l_evaluations = compute_pseudorandom_linear_combination(mtree, G2, steps, precision, d_evaluations, p_evaluations, b_evaluations)
  l_mtree = merkelize(l_evaluations)

  branches = compute_merkle_spot_checks(mtree, l_mtree, precision, skips)

  # Return the Merkle roots of P and D, the spot check Merkle
  # proofs, and low-degree proofs of P and D
  o = [
      mtree[1], l_mtree[1], branches,
      prove_low_degree(
          l_evaluations,
          G2,
          # TODO(rbharath): Why is this 2x?
          steps * constraint_degree,
          modulus,
          exclude_multiples_of=extension_factor)
  ]
  print("STARK computed in %.4f sec" % (time.time() - start_time))
  return o

def verify_proof_at_position(inp, output, ks, G2, steps, skips, skips2, precision, proof, i, pos, last_step_position, constants_polynomials, step_fn):
  """Verifies merkle proof at given position in extended trace"""
  k1, k2, k3, k4 = ks
  m_root, l_root, branches, fri_proof = proof
  x = f.exp(G2, pos)
  # TODO(rbharath): Why does exponentiating to steps make
  # sense here?  I think this is to compute the pseudorandom
  # linear combination.
  x_to_the_steps = f.exp(x, steps)
  # TODO(rbharath): Why do i*3, i*3+1, i*3+2 make sense?
  mbranch1 = verify_branch(m_root, pos, branches[i * 3])
  mbranch2 = verify_branch(m_root, (pos + skips) % precision,
                            branches[i * 3 + 1])
  l_of_x = verify_branch(l_root, pos, branches[i * 3 + 2], output_as_int=True)

  p_of_x = int.from_bytes(mbranch1[:32], 'big')
  p_of_g1x = int.from_bytes(mbranch2[:32], 'big')
  d_of_x = int.from_bytes(mbranch1[32:64], 'big')
  b_of_x = int.from_bytes(mbranch1[64:], 'big')

  zvalue = f.div(f.exp(x, steps) - 1, x - last_step_position)
  k_of_xs = []
  for constants_mini_polynomial in constants_polynomials:
    k_of_x = f.eval_poly_at(constants_mini_polynomial, f.exp(x, skips2))
    k_of_xs.append(k_of_x)

  # Check transition constraints C(P(x)) = Z(x) * D(x)
  assert (p_of_g1x - step_fn(f, p_of_x, k_of_xs) - zvalue * d_of_x) % modulus == 0

  # Check boundary constraints B(x) * Q(x) + I(x) = P(x)
  interpolant = f.lagrange_interp_2([1, last_step_position], [inp, output])
  zeropoly2 = f.mul_polys([-1, 1], [-last_step_position, 1])
  assert (p_of_x - b_of_x * f.eval_poly_at(zeropoly2, x) - f.eval_poly_at(
      interpolant, x)) % modulus == 0

  # Check correctness of the linear combination
  assert (l_of_x - d_of_x - k1 * p_of_x - k2 * p_of_x * x_to_the_steps -
          k3 * b_of_x - k4 * b_of_x * x_to_the_steps) % modulus == 0


def verify_proof(inp, steps, constants, output, proof, step_fn,
    constraint_degree=2):
  """Verifies a STARK
  
  Parameters
  ----------
  constraint_degree: int
    The degree of the constraint being considered
  """
  m_root, l_root, branches, fri_proof = proof
  start_time = time.time()
  assert steps <= 2**32 // extension_factor
  # ALl constants should be of same length so we check the first
  #assert is_a_power_of_2(steps) and is_a_power_of_2(len(constants[0]))
  assert is_a_power_of_2(steps)
  assert len(constants[0]) <= steps

  precision = steps * extension_factor

  # Get (steps)th root of unity
  G2 = f.exp(7, (modulus - 1) // precision)
  skips = precision // steps

  # Gets the polynomial representing the constants
  skips2 = steps // len(constants[0])

  deg = len(constants)
  constants_polynomials = []
  constants_extensions = []
  for d in range(deg):
    deg_constants = constants[d]
    skips2 = steps // len(deg_constants)
    constants_mini_polynomial = fft(
        deg_constants, modulus, f.exp(G2, extension_factor * skips2), inv=True)
    constants_mini_extension = fft(constants_mini_polynomial, modulus,
                                  f.exp(G2, skips2))
    assert len(constants_mini_extension) == precision // skips2
    constants_polynomials.append(constants_mini_polynomial)
    constants_extensions.append(constants_mini_extension)

  # Verifies the low-degree proofs
  assert verify_low_degree_proof(
      l_root,
      G2,
      fri_proof,
      steps * constraint_degree,
      modulus,
      exclude_multiples_of=extension_factor)

  # Performs the spot checks
  k1 = int.from_bytes(blake(m_root + b'\x01'), 'big')
  k2 = int.from_bytes(blake(m_root + b'\x02'), 'big')
  k3 = int.from_bytes(blake(m_root + b'\x03'), 'big')
  k4 = int.from_bytes(blake(m_root + b'\x04'), 'big')
  ks = [k1, k2, k3, k4]
  samples = spot_check_security_factor
  positions = get_pseudorandom_indices(
      l_root, precision, samples, exclude_multiples_of=extension_factor)
  last_step_position = f.exp(G2, (steps - 1) * skips)
  for i, pos in enumerate(positions):
    verify_proof_at_position(inp, output, ks, G2, steps, skips, skips2, precision, proof, i, pos, last_step_position, constants_polynomials, step_fn)

  print('Verified %d consistency checks' % spot_check_security_factor)
  print('Verified STARK in %.4f sec' % (time.time() - start_time))
  return True
