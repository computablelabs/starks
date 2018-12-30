"""This file holds the classes necessary to implement the algebraic placement
and routing (APR) reduction. This reduction is a crucial step on the way to
zero-knowledge support since the prover applies randomness during the APR
reduction to achieve a perturbed version of the original AIR representation.
"""

class APR(object):
  """A class holding an instance of the APR problem.

  Recall that instance of the APR problem is a tuple

    x = (F, T, N, Phi, L, Lcmp, rhovec, rhocomp)
    
  Let's define each of these terms in sequence.

  - F is a finite field of characteristic 2 (TODO(rbharath): Can this be changed to characteristic p?)
  - T is a set of indices (TODO(rbharath): what are the indices?)
  - N is a subset of T x Aff_1(F). Recall Aff_1(F) is the set of linear polynomials over F.
  - Phi is  subset of (F x F^N)->F. That is, it's a set of functions over the
    variables {X_loc, X_n}_{n in N}.
  - L is an F_2 affine subspace of F (TODO(rbharath): represented how?)
  - Lcmp is another F_2 affine subspace of F (TODO(rbharath): represented how?)
  - phovec in (0, 1)^T is a list of rates between 0 and 1
  - phocmp in (0, 1) is a rate between 0 and 1.
  """
  def __init__(self):
    pass

  def get_witness(self):
    """A witness w^hat is in (L^F)^T. That is, it's a set of functions indexed
    by set T. Each function maps from space L to field F."""
    return None
