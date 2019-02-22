"""This file holds the classes necessary to implement the algebraic placement
and routing (APR) reduction. This reduction is a crucial step on the way to
zero-knowledge support since the prover applies randomness during the APR
reduction to achieve a perturbed version of the original AIR representation.

This file transforms the trace into a path (?) in an affine graph. Here, an affine graph has elements of finite field F (recall these correspond to "words" in a processor architecture) and the affine graph is a directed graph with vertices F and specified edges.
"""

import math
from starks.air import AIR 
from starks.reedsolomon import AffineSpace
from starks.polynomial import polynomials_over
from starks.poly_utils import generate_primitive_polynomial
from starks.multivariate_polynomial import multivariates_over
from starks.poly_utils import construct_affine_vanishing_polynomial
from starks.poly_utils import lagrange_interp
from starks.poly_utils import draw_random_interpolant

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

  TODO(rbharath): This class has a lot of really gnarly methods that should be
  factored out into more general functions to minimize the amount of state
  being handled.

  TODO(rbharath): The low-degree extension (additive FFT) is used in this
  class! Where are the specific incantions?
  """
  def __init__(self, air: AIR, zero_knowledge_expansion=5):
    """Initiates the APR object using an AIR object.

    Implements the AIR->APR transform from the STARKs paper.
    """
    self.air = air 
    modulus = air.field.p
    self.width = air.width
    # field = (Z/2[g]/h(g))
    self.field = air.field
    # base_field = Z/2
    base_field = self.field.base_field
    self.basePolys = polynomials_over(base_field)
    # h(g)
    h = self.field.ideal_generator
    # g
    g = self.basePolys([0, 1])

    # Some constants
    T = air.steps
    self.t = int(math.log(T, 2))
    # chosen so deg(C) <= 2^d
    # Setting to arbitrary value for now.
    # TODO(rbharath): What is the right value of this?
    self.d = 10
    # TODO(rbharath): How should this constant be set correctly?
    self.R = 5
    # This is the zero-knowledge expansion
    self.k = zero_knowledge_expansion

    self.polysOver = polynomials_over(self.field).factory
    self.Tau = list(range(self.width))
    # Element of (Z/2[g]/h(g))
    self.zeta = generate_primitive_polynomial(modulus, self.width)
    # Neighbors
    self.Nbrs = self.construct_neighbors(self.Tau, self.zeta, g, self.polysOver)

    # Define the affine spaces (TODO(rbharath): Fill out stub)
    self.H = AffineSpace(base_field, [g**k for k in range(self.t)])
    self.H0 = AffineSpace(base_field, [g**k for k in range(self.t-1)])
    self.H1 = AffineSpace(base_field, [g**k for k in range(self.t-1)], g**(self.t-1))
    self.L = self.construct_L(g) 
    self.Lcmp = self.construct_L_cmp(g) 
    self.Z_boundaries = self.construct_Z_boundaries(air.B)
    self.Eps_boundaries = self.construct_Eps_boundaries(air.B)
    self.rho_js = self.compute_rho_js(self.Z_boundaries, self.L)
    self.rho_cmp = self.compute_rho_cmp(self.Lcmp)

    # X_loc + {X_N}_{n in Nbrs}
    num_Phi_vars = 1 + len(self.Nbrs)
    PhiPolys = multivariates_over(self.field, num_Phi_vars).factory
    self.Phi = self.construct_Phi_polynomials(air, PhiPolys, g, self.zeta)

  def tilde_expansion(indices, neighbor):
    """Performs the Tilde expansion of a neighbor.

    The STARKs paper defines a "tilde" expansion as follows:

    Tilde{(tau, N)} = X_{(tau, N)} * Z_{B, tau}(N(X_loc)) + E_{B, tau}(N(X_loc))

    This has an expansion to sequences as follows

    Tilde{(tau_1,...,tau_n, N)} = Tilde{(tau_1, N)},...,Tilde{(tau_n, N)}
    """
    # TODO(rbharath): Fill this out.
    return []

  def construct_L(self, g):
    """Constructs the affine space L"""
    basis = [g**i for i in range(1 + self.k + self.R + self.t + self.d)] + [g**(1+self.k+self.R+self.t+self.d)*(1+g)]
    shift = g**(1+self.k+self.R+self.t+self.d)
    return AffineSpace(self.field, basis, shift)

  def construct_L_cmp(self,g):
    """Constructs the affine space Lcmp"""
    basis = [g**i for i in range(1 + self.k + self.R + self.t + self.d)]
    shift = g**(1+self.k+self.R+self.t+self.d)
    return AffineSpace(self.field, basis, shift)

  def compute_rho_js(self, Z_boundaries, L):
    """Computes rho_j

    rho_j = (2**(k+t) - deg(Z_{B_j}))/|L|
    """
    rho_js = []
    for Z_B_j in Z_boundaries:
      # TODO(rbharath): Using len(L) converts to C-ints somewhere which causes overflow problems
      rho_j = int(2**(self.k+self.t) - Z_B_j.degree())/L.__len__()
      rho_js.append(rho_j)
    return rho_js

  def compute_rho_cmp(self, Lcmp):
    """Computes rho_cmp

    rho_cmp = (1 + 2**(k+t+d))/|Lcmp|
    """
    # TODO(rbharath): Using len(L) converts to C-ints somewhere which causes overflow problems
    return (1 + 2**(self.k+self.t+self.d))/Lcmp.__len__()


  def construct_Phi_polynomials(self, comp, PhiPolys, g, zeta):
    """Constructs Phi polynomials that encode transition and boundary constraints.

    For each P in comp.Polys:
      Phi_P_0 = (X_{loc}*(X_{loc} - 1))/Z_H0 * P(Tilde{(1,...,w, n^{id})}, Tilde{(1,...,w, n_0^cyc)})
      Phi_P_1 = 1/Z_H0 * P(Tilde{(1,...,w, n^{id})}, Tilde{(1,...,w, n_1^cyc)})
    """
    Phis = []
    Z_H0 = construct_affine_vanishing_polynomial(self.field, self.H0)
    # Define X_loc
    X_loc = PhiPolys({(1,) + (0,)*len(self.Nbrs): 1})
    # n_id(x) = x
    n_id = self.polysOver([0, 1])
    # n_cyc_1 = gx + zeta
    n_cyc_1 = self.polysOver([zeta, g])
    # n_cyc_0 = gx
    n_cyc_0 = self.polysOver([0, g])
    for P in comp.Polys:
      # TODO(rbharath): How does this division work? Division in multivariate
      # polynomial rings is gnarly to get right.
      # TODO(rbharath): Does evaluation on other polynomials work out of the box?
      TODO = 1
      Phi_P_0 = ((X_loc * (X_loc - 1))/Z_H0) * P(
        self.tilde_expansion(range(self.width), n_id),
        self.tilde_expansion(range(self.width), n_0_cyc))
      Phi_P_1 = 1/Z_H0 * P(
        self.tilde_expansion(range(self.width), n_id),
        self.tilde_expansion(range(self.width), n_1_cyc))
      Phis.append(Phi_P_0)
      Phis.append(Phi_P_1)
    return Phis


  # TODO(rbharath): This is broken. Fix!!
  def construct_Z_boundaries(self, B):
    """Constructs Z_{B,j}(x) boundary constraint polynomials

    The Z_{B,j} polynomial is defined in the STARKs paper by the following equation

    Z_{B,j}(x) = \prod_{(i, j, alpha) \in B}
    """
    ##########################################
    print("B")
    print(B)
    ##########################################
    accums = []
    x = self.polysOver([0, 1])
    g = self.basePolys([0, 1])
    ##########################################
    print("self.width")
    print(self.width)
    print("self.zeta")
    print(self.zeta)
    print("(g**2) % self.zeta")
    print((g**2) % self.zeta)
    ##########################################
    for w in range(self.width):
      accum = self.polysOver([self.basePolys(1)])
      for (i, j, alpha) in B:
        # x - g**i % zeta
        term = x - self.polysOver([((g**i) % self.zeta)])
        accum = accum * term
      accums.append(accum)
    return accums

  def construct_Eps_boundaries(self, B):
    """Construcs E_{B,j}(x) boundary constraint polynomials.

    The E_{B,j} polynomial is defined in the STARKs paper by the following equation.

    E_{B,j}(x) = interpolant([g^i % zeta: alpha for (i, j, alpha) in B])

    TODO(rbharath): This function doesn't work since lagrange interp doesn't
    currently work over finite fields.
    """
    accums = []
    g = self.basePolys([0, 1])
    for w in range(self.width):
      xs = [(g**i) % self.zeta for (i, _, _) in B]
      ys = [alpha for (_, _, alpha) in B]
      # TODO(rbharath): This interpolation call is broken!
      interp = lagrange_interp(self.field, xs, ys)
      accums.append(interp)
    return accums

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

  def get_boundary_conditions(self):
    """Retrieves boundary conditions B."""
    raise NotImplementedError

  def generate_witness(self):
    """A witness w^hat is in (L^F)^T. That is, it's a set of functions indexed
    by set T. Each function maps from space L to field F.
   
    TODO(rbharath): I think the reason zero-knowledge holds here is that it's
    hard to find g**i % zeta for the verifier since the only thing they're
    given is the coefficients of a high dimensional interpolant polynomial.
    
    """
    comp = self.comp
    Z_boundary = self.Z_boundaries
    E_boundary = self.Eps_boundaries
    w_air = comp.get_witness()
    w_apr = []
    g = self.basePolys([0, 1])
    for (w_air_j, Z_B_j, E_B_j) in zip(w_air, Z_boundary, E_boundary):
      xs = [g**i % self.zeta for i in range(2**self.t)]
      ys = [w_air_j[i] for i in range(2**self.t)]
      Q_j = draw_random_interpolant(2**(self.k+self.t), xs, ys)
      # TODO(rbharath): How does this division work? These are multivariate polynomials I think.
      w_apr_j = (Q_j - E_B_j)/Z_B_j
      w_apr.append(w_apr_j)
    return w_apr
