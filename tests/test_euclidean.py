import pytest
#from starks.euclidean import gcd
#from starks.euclidean import extended_euclidean_algorithm
#from starks.modp import IntegersModP


#class TestEuclidean(unittest.TestCase):
#  """"
#  Basic tests for the extended euclidean algorithm
#  """
#
#  def test_gcd_ints(self):
#    """Test that the GCD functions correctly on ints."""
#    assert gcd(7, 9) == 1
#    assert gcd(8, 18) == 2
#    assert gcd(-12, 24) == -12
#    # gcd is only unique up to multiplication by a unit, and so sometimes we'll get negatives.
#    assert gcd(12, -24) == 12
#    assert gcd(4864, 3458) == 38
#
#  def test_extended_euclidean_ints(self):
#    """Test that the extended Euclidean algorithm works correctly on ints."""
#    assert extended_euclidean_algorithm(4864, 3458) == (32, -45, 38)
#    assert 32 * 4864- 45 * 3458 == 38
#    assert extended_euclidean_algorithm(3458, 4864) == (-45, 32, 38)
#
#  def test_mod_2_gcd(self):
#    """Test that the GCD works in mod-2 arithmetic."""
#    Mod2 = IntegersModP(2)
#    assert Mod2(1) == gcd(Mod2(1), Mod2(0))
#    assert Mod2(1) == gcd(Mod2(1), Mod2(1))
#    assert Mod2(0) == gcd(Mod2(2), Mod2(2))
#
#  def test_mod_7_gcd(self):
#    """Test that the GCD works in mod-7 arithmetic."""
#    Mod7 = IntegersModP(7)
#    # TODO(rbharath): Why are these modular equations right? Is there a way to
#    # do simple mental arithmetic to calculate these values?
#    assert Mod7(6) == gcd(Mod7(6), Mod7(14))
#    assert Mod7(2) == gcd(Mod7(6), Mod7(9))
#
#  def test_mod_large_gcd(self):
#    """Test that GCD works with a larger prime."""
#    ModHuge = IntegersModP(9923)
#    assert ModHuge(38) == gcd(ModHuge(4864), ModHuge(3458))
#    assert (ModHuge(32), ModHuge(-45), ModHuge(38)) == extended_euclidean_algorithm(
#      ModHuge(4864), ModHuge(3458))
#
