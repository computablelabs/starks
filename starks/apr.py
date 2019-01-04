"""This file holds the classes necessary to implement the algebraic placement
and routing (APR) reduction. This reduction is a crucial step on the way to
zero-knowledge support since the prover applies randomness during the APR
reduction to achieve a perturbed version of the original AIR representation.
"""

import math
from starks.air import Computation
from starks.polynomial import polynomials_over
from starks.poly_utils import generate_primitive_polynomial

class APR(object):
  """A class holding an instance of the APR problem.

  Recall that instance of the APR problem is a tuple

    x = (F, Tau, N, Phi, L, Lcmp, rhovec, rhocomp)
    
  Let's define each of these terms in sequence.

  - F is a finite field of characteristic 2 (TODO(rbharath): Can this be changed to characteristic p?)
  - Tau is a set of indices, typically [w]
  - N is a subset of Tau x Aff_1(F). Recall Aff_1(F) is the set of linear polynomials over F.
  - Phi is  subset of (F x F^N)->F. That is, it's a set of functions over the
    variables {X_loc, X_n}_{n in N}.
  - L is an F_2 affine subspace of F (TODO(rbharath): represented how?)
  - Lcmp is another F_2 affine subspace of F (TODO(rbharath): represented how?)
  - phovec in (0, 1)^Tau is a list of rates between 0 and 1
  - phocmp in (0, 1) is a rate between 0 and 1.
  """
  def __init__(self, comp: Computation):
    """Initiates the APR object using an AIR Computation object.

    Implements the AIR->APR transform from the STARKs paper.
    """
    modulus = comp.field.p
    width = comp.dims
    # field = (Z/2[g]/h(g))
    field = comp.field
    # base_field = Z/2
    base_field = field.field
    basePolys = polynomials_over(base_field)
    # h(g)
    h = field.ideal_generator
    # g
    g = basePolys([0, 1])

    # Some constants
    self.t = math.log(T, 2)
    # chosen so deg(C) <= 2^d
    # Setting to arbitrary value for now.
    # TODO(rbharath): What is the right value of this?
    self.d = 10

    polysOver = polynomials_over(field).factory
    self.Tau = list(range(width))
    # Element of (Z/2[g]/h(g))
    zeta = generate_primitive_polynomial(modulus, width)
    # Neighbors
    N = self.construct_neighbors(self.Tau, zeta, g, polysOver)

  def construct_z_B_j(self):
    """Constructs Z_{B,j}(x) boundary constraint polynomial

    The Z_{B,j} polynomial is defined in the STARKs paper by the following equation

    Z_{B,j}(x) = \prod_{(i, j, alpha) \in B}
    """
    pass

  def construct_neighbors(self, Tau, zeta, g, polysOver):
    """Helper method to construct neighbor set."""
    neighbors = []
    for tau in Tau:
      # n_id(x) = x
      n_id = polysOver([0, 1])
      # n_cyc_1 = gx + zeta
      n_cyc_1 = polysOver([zeta, g])
      # n_cyc_0 = gx
      n_cyc_0 = polysOver([0, g])
      neighbors.extend([(tau, n_id), (tau, n_cyc_1), (tau, n_cyc_0)])
    return neighbors


  def generate_witness(self):
    """A witness w^hat is in (L^F)^T. That is, it's a set of functions indexed
    by set T. Each function maps from space L to field F."""
    return None
