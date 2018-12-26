import unittest
from starks.modp import IntegersModP
from starks.rationals_modp import RationalsModP

class TestModP(unittest.TestCase):
  """Basic tests for Mod-p numbers."""

  def test_construction(self):
    """Test that you can initialize a mod-p number."""
    mod7 = IntegersModP(7)

  def test_arithmetic(self):
    """Test Basic arithmetic."""
    mod7 = IntegersModP(7)
    # Check basic equality
    assert mod7(5) == mod7(5)
    assert mod7(5) == mod7(12)

    # Note this works since the @typecheck decorator silently casts the 1 to an IntegerModP 
    assert mod7(5) == 1 / mod7(3) # 3 * 5 = 15 == 1 mod 7
    assert mod7(1) == mod7(3) * mod7(5)
    # Check that the @typecheck decorator works correctly
    assert mod7(3) == mod7(3) * 1
    assert mod7(2) == mod7(5) + mod7(4)

    assert mod7(0) == mod7(3) + mod7(4)

  def test_large_modulus(self):
    """Runs basic tests in large modulus needed for starks."""
    modulus = 2**256 - 2**32 * 351 + 1
    modM = IntegersModP(modulus)
    assert modM(5) != modM(11)
    assert modM(2**32) != modM(2**64)
    assert modM(2**64) != modM(2**128)
    assert modM(2**256) == modM(2**32 * 351 - 1)

  def test_rationals_construction(self):
    """Test that rationals can be constructed."""
    modulus = 3
    ratmod3 = RationalsModP(modulus)
    rational = ratmod3(1, 2) # 1/2

  def test_rationals_arithmetic(self):
    """Test that rational arithmetic is sensible."""
    modulus = 31 
    ratmod = RationalsModP(modulus)
    assert ratmod(1, 2) * ratmod(1, 2) == ratmod(1, 4) 
    assert ratmod(3, 4) * ratmod(5, 6) == ratmod(5, 8) 
    assert 1 / ratmod(3, 4) == ratmod(4, 3) 
    assert 1 + ratmod(1, 2) == ratmod(3, 2)
    assert -ratmod(1, 2) == ratmod(-1, 2)
    assert 1 - ratmod(1, 2) == ratmod(1, 2)

  def test_count_unique_rationals(self):
    """Test count of unique number of rational numbers
    
    |Q/p| == p by some abstract algebra. Check this maps out in our implementation. 
    """
    modulus = 7
    ratmod = RationalsModP(modulus)
    uniques = []
    for num in range(modulus):
      for den in range(1, modulus):
        cur = ratmod(num, den)
        is_unique = True 
        for unique in uniques:
          if cur == unique:
            is_unique = False
        if is_unique:
          uniques.append(cur)
    assert len(uniques) == modulus

