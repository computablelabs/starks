"""This file holds the classes needed to define the algebraic linking IOP (ALI) protocol.

"""

from starks.reedsolomon import ReedSolomonCode

class ALIProver(object):
  """The ZK-stark prover for the ALI protocol.
  
  TODO(rbharath): At present the protocol implemented is interactive. Implement
  the Fiat-shamir heuristic so we get a single ZK-stark proof. 
  """

  def __init__(self, apr):
    """TODO(rbharath): Fill this out."""
    self.apr = apr
    self.rs_f = ReedSolomonCode(apr.field, apr.L, apr.rhomax)
    self.rs_g = ReedSolomonCode(apr.field, apr.Lcmp, apr.rhocmp)

  def generate_proof(self):
    """Generates ZK-proof."""
    f_mask = self.rs_f.draw_random_sample()
    g_mask = self.rs_g.draw_random_sample()
    w_apr = self.apr.generate_witness()
    return (w_apr, f_mask, g_mask)

  def generate_random_constraint(self, rphi):
    """Generate the random constraint.

    phi_R(alpha) = \sum_{phi in Phi} r_phi * phi(alpha)
    """
    Phi = self.apr.Phi
    phi_R = 0
    for i, phi in enumerate(Phi):
      phi_R += r[i] * phi
    return phi_R


class APIVerifier(object):
  """The verifier for the ZK-stark protocol."""
  pass
  

