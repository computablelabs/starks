"""This file holds classes necessary to implement Reed-Solomon codes"""
import itertools
from typing import List
from starks.numbertype import Field
from starks.numbertype import FieldElement
from starks.numbertype import Poly

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
      print("basis_vals")
      print(basis_vals)
      elt = self.shift
      for val, basis_elt in zip(basis_vals, self.basis):
        elt += val * basis_elt
      print("elt")
      print(elt)
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
