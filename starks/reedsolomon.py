"""This file holds classes necessary to implement Reed-Solomon codes"""
from typing import List
from starks.numbertype import Field
from starks.numbertype import FieldElement

class SmoothMultiplicativeGroup(object):
  """Defines a smooth multiplicative group.

  A smooth multiplicative group is a group of order 2^n which
  is a multiplicative subgroup of some finite field.
  """

  def __init__(self):
    pass

class AffineSpace(object):
  """Defines an affine space (of polynomials typically)."""
  def __init__(self, field: Field, basis: List[FieldElement], shift=None):
    self.field = field
    self.basis = basis
    if shift is None:
      self.shift = self.field(0)
    else:
      self.shift = shift

  def __len__(self):
    """Returns the size of the affine space."""
    # TODO(rbharath): Implement this in a reasonable fashion
    field_size = int(self.field.p)**self.field.m
    basis_len = int(len(self.basis))
    length = field_size**basis_len
    return length

  def __iter__(self):
    """Iterates over the elements of the affine space."""
    field_iterators = []
    for _ in range(len(self.basis)):
      field_iterators.append(self.field.__iter__())
    for basis_vals in itertools.product(*field_iterators):
      elt = self.shift
      for val, basis_elt in zip(basis_vals, self.basis):
        elt += val * basis_elt
      yield elt


class ReedSolomonCode(object):
  """Defines a Reed Solomon Code.
  
  A Reed Solomon code RS[F, S, rho] is defined by the following parameters

  - F: A finite field
  - S: A subset of the finite field
  - rho: In (0, 1) is a rate parameter

  RS[F, S, rho] is the family of functions f:S -> F that are evaluations of
  poynomials of degree < rho |S|.

  Note also that the triples x_{RS} = (F, S, rho) are said to be instances of the
  RS Proximity testing problem (RPT). We say w_{RS} is a witness for x_{RS} if
  it is a function w_{RS}: S -> F. We say that w_{RS} satisfies x_{RS} if and
  only if w_{RS} \in RS[F, RS, rho]. The relation R_{RPT} is the set of pair
  (x_{RS}, w_{RS}).

  There is a special case in which"

  - F is a binary field
  - S is an affine coset of a F_2 linear subspace of F. That is, S = {a_0 + \sum_{i=1}^k alpha_i a_i} where (a_1,...,a_k) is a set of k-linearly independent elements.
  - rho = 2^{-Rcurly} where Rcurly is a positive integer.

  If these conditions are met, we saw that this defines the binary RPT relation
  R_{BRPT}. In this case, RS[F, S, rho] is also called a binary additive RS
  code family.

  In this class, we assume that only binary RPT problems are defined.
  """

  def __init__(self, field: Field, S: List[FieldElement], rho: float):
    """Initializes the Reed Solomon Code."""
    self.field = field
    self.S = S
    self.rho = rho

  def draw_random_sample(self):
    """Draws a random sample from this RS code."""
    raise NotImplementedError


class FRI(object):
  """Implements Fast Reed Solomon Interactive Oracle Protocol"""
  def __init__(self, rs: ReedSolomon):
    self.rs = rs

  def prove_proximity(f: Poly,
                      root_of_unity: FieldElement,
                      maxdeg_plus_1: int,
                      exclude_multiples_of:int = 0,
                      fri_spot_check_security_factor:int = 40):
    """
    Generate an FRI proof that the polynomial that has the
    specified values at successive powers of the specified root
    of unity has a degree lower than maxdeg_plus_1
    
    We use maxdeg+1 instead of maxdeg because it's more
    mathematically convenient in this case.

    Note that if values is a n-degree polynomial, root_of_unity
    should be a n-th root of unity.
    """
    # Is this the right iteration?
    values = [f(x) for x in rs]
    # If the degree we are checking for is less than or equal
    # to 32, use the polynomial directly as a proof
    # TODO(rbharath): Why does this make sense?
    if maxdeg_plus_1 <= 16:
      print('Produced FRI proof')
      return [[x.to_bytes() for x in values]]

    # Calculate the set of x coordinates
    xs = get_power_cycle(root_of_unity, field)
    
    assert len(values) == len(xs)

    # Put the values into a Merkle tree. This is the root that
    # the proof will be checked against. Note this is a list of
    # length 2*len(values) storing the complete merkle tree.
    m = merkelize(values)

    # Select a pseudo-random x coordinate
    # This is the merkle-root of the polynomial.
    #special_x = int.from_bytes(m[1], 'big') % modulus
    special_x = field(m[1])

    # Calculate the "column" at that x coordinate (see
    # https://vitalik.ca/general/2017/11/22/starks_part_2.html)
    # We calculate the column by Lagrange-interpolating each
    # row, and not directly from the polynomial, as this is more
    # efficient
    quarter_len = len(xs) // 4
    x_polys = multi_interp_4(field,
        [[xs[i + quarter_len * j] for j in range(4)] for i in range(quarter_len)],
        [[values[i + quarter_len * j]
          for j in range(4)]
        for i in range(quarter_len)])
    #column = [f.eval_quartic(p, special_x) for p in x_polys]
    column = [p(special_x) for p in x_polys]
    m2 = merkelize(column)

    # Pseudo-randomly select y indices to sample
    ys = get_pseudorandom_indices(
        m2[1], len(column), fri_spot_check_security_factor, exclude_multiples_of=exclude_multiples_of)

    # Compute the Merkle branches for the values in the
    # polynomial and the column
    branches = []
    for y in ys:
      branches.append([mk_branch(m2, y)] +
                      [mk_branch(m, y + (len(xs) // 4) * j) for j in range(4)])

    o = [m2[1], branches]

    # Recurse...
    return [o] + self.prove_proximity(
        column,
        root_of_unity**4,
        maxdeg_plus_1 // 4,
        field,
        exclude_multiples_of=exclude_multiples_of)

  def verify_proximity(self, f) -> bool:
    """Verifies proximity of this function this RS code."""
    raise NotImplementedError
