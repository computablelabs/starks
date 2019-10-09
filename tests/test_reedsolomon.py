import pytest
#from starks.polynomial import polynomials_over
#from starks.modp import IntegersModP
#from starks.finitefield import FiniteField
#from starks.reedsolomon import AffineSpace

#class TestAffineSpace(unittest.TestCase):
#  """
#  Basic tests for Affine Space classes 
#  """
#
#  def test_affine_space_construction(self):
#    """Construct affine space."""
#    p = 2
#    m = 1
#    # Degree of space we construct
#    t = 2
#    Zp = IntegersModP(p)
#    basePolys = polynomials_over(Zp)
#    # g
#    g = basePolys([0, 1])
#    field = FiniteField(p, m)
#    H0 = AffineSpace(Zp, [g**k for k in range(t)])
#
#  def test_affine_space_iteration(self):
#    p = 2
#    m = 2
#    # Degree of space we construct
#    t = 2
#    Zp = IntegersModP(p)
#    basePolys = polynomials_over(Zp)
#    # g
#    g = basePolys([0, 1])
#    field = FiniteField(p, m)
#    H0 = AffineSpace(Zp, [g**k for k in range(t)])
#    elts = [aff for aff in H0]
#    #########################################
#    print("elts")
#    print(elts)
#    #assert len(elts) == (p)**t
#    assert len(elts) == 4 
#    #########################################
