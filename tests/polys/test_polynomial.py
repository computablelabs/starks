import pytest
from starks.finitefield import IntegersModP
from fractions import Fraction
from starks.polynomial import polynomials_over

def test_equality():
    """Basic test of polynomial equality."""
    Mod5 = IntegersModP(5)
    Mod11 = IntegersModP(11)
  
    polysOverQ = polynomials_over(Fraction).factory
    polysMod5 = polynomials_over(Mod5).factory
    polysMod11 = polynomials_over(Mod11).factory
    for p in [polysOverQ, polysMod5, polysMod11]:
        # equality
        assert p([]) == p([])
        assert p([1,2]) == p([1,2])
        assert p([1,2,0]) == p([1,2,0,0])

#def test_addition():
#  """Basic test of polynomial addition."""
#  Mod5 = IntegersModP(5)
#  Mod11 = IntegersModP(11)
#
#  polysOverQ = polynomials_over(Fraction).factory
#  polysMod5 = polynomials_over(Mod5).factory
#  polysMod11 = polynomials_over(Mod11).factory
#  for p in [polysOverQ, polysMod5, polysMod11]:
#   # addition
#   assert p([1,2,3]) == p([1,0,3]) + p([0,2])
#   assert p([1,2,3]) == p([1,2,3]) + p([])
#   assert p([5,2,3]) == p([4]) + p([1,2,3])
#   assert p([1,2]) == p([1,2,3]) + p([0,0,-3])
#
#def test_subtraction():
#  """Basic test of polynomial subraction."""
#  Mod5 = IntegersModP(5)
#  Mod11 = IntegersModP(11)
#
#  polysOverQ = polynomials_over(Fraction).factory
#  polysMod5 = polynomials_over(Mod5).factory
#  polysMod11 = polynomials_over(Mod11).factory
#  for p in [polysOverQ, polysMod5, polysMod11]:
#   # subtraction
#   assert p([1,-2,3]) == p([1,0,3]) - p([0,2])
#   assert p([1,2,3]) == p([1,2,3]) - p([])
#   assert p([-1,-2,-3]) == p([]) - p([1,2,3])
#
#def test_multiplication():
#  """Basic test of polynomial multiplication."""
#  Mod5 = IntegersModP(5)
#  Mod11 = IntegersModP(11)
#
#  polysOverQ = polynomials_over(Fraction).factory
#  polysMod5 = polynomials_over(Mod5).factory
#  polysMod11 = polynomials_over(Mod11).factory
#  for p in [polysOverQ, polysMod5, polysMod11]:
#   # multiplication
#   assert p([1,2,1]) == p([1,1]) * p([1,1])
#   assert p([2,5,5,3]) == p([2,3]) * p([1,1,1])
#   assert p([0,7,49]) == p([0,1,7]) * p([7])
#
#def test_division():
#  """Basic test of polynomial division."""
#  Mod5 = IntegersModP(5)
#  Mod11 = IntegersModP(11)
#
#  polysOverQ = polynomials_over(Fraction).factory
#  polysMod5 = polynomials_over(Mod5).factory
#  polysMod11 = polynomials_over(Mod11).factory
#  for p in [polysOverQ, polysMod5, polysMod11]:
#    # division
#    assert p([1,1,1,1,1,1]) == p([-1,0,0,0,0,0,1]) / p([-1,1])
#    assert p([-1,1,-1,1,-1,1]) == p([1,0,0,0,0,0,1]) / p([1,1])
#    assert p([]) == p([]) / p([1,1])
#    assert p([1,1]) == p([1,1]) / p([1])
#    assert p([1,1]) == p([2,2]) / p([2])
#
#def test_modulus():
#  """Basic test of polynomial modulus."""
#  Mod5 = IntegersModP(5)
#  Mod11 = IntegersModP(11)
#
#  polysOverQ = polynomials_over(Fraction).factory
#  polysMod5 = polynomials_over(Mod5).factory
#  polysMod11 = polynomials_over(Mod11).factory
#  for p in [polysOverQ, polysMod5, polysMod11]:
#    # modulus
#    assert p([]) == p([1,7,49]) % p([7])
#    assert p([-7]) == p([-3,10,-5,3]) % p([1,3])
#
#def test_division_more():
#  """More division tests"""
#  Mod5 = IntegersModP(5)
#  Mod11 = IntegersModP(11)
#
#  polysOverQ = polynomials_over(Fraction).factory
#  polysMod5 = polynomials_over(Mod5).factory
#  polysMod11 = polynomials_over(Mod11).factory
#
#  assert polysOverQ([Fraction(1,7), 1, 7]) == polysOverQ([1,7,49]) / polysOverQ([7])
#  assert polysMod5([1 / Mod5(7), 1, 7]) == polysMod5([1,7,49]) / polysMod5([7])
#  assert polysMod11([1 / Mod11(7), 1, 7]) == polysMod11([1,7,49]) / polysMod11([7])
#
#def test_polynomial_call():
#  """Test evaluation of polynomials."""
#  mod5 = IntegersModP(5)
#  polysMod5 = polynomials_over(mod5).factory
#  # 1 + x
#  poly = polysMod5([1, 1])
#  # z = 3
#  z = mod5(3)
#  assert z + 1 == poly(z)
#  # 1 + x + x^2 (1 + 3 + 9 == 13 == 3)
#  poly2 = polysMod5([1, 1, 1])
#  assert 1 + z + z**2 == poly2(z)
#  assert poly2(z) == mod5(3)
#
#
