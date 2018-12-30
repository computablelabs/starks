import unittest
from starks.modp import IntegersModP
from starks.poly_utils import zpoly
from starks.poly_utils import multi_inv
from starks.poly_utils import lagrange_interp
from starks.poly_utils import lagrange_interp_2
from starks.poly_utils import lagrange_interp_4
from starks.poly_utils import multi_interp_4 
from starks.polynomial import polynomials_over
from starks.utils import get_power_cycle

class TestPolyUtils(unittest.TestCase):
  """
  Basic tests for polynomial utility functions. 
  """

  def test_basic(self):
    """Basic test"""
    #field7 = PrimeField(7)
    mod7 = IntegersModP(7)
    # 12 % 7 == 5
    assert mod7(6) + mod7(6) == mod7(5)

    # 6^-1 = 6
    assert 1/mod7(6) == mod7(6)

  def test_zpoly(self):
    """Test construction of polynomials with specified root"""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysMod7 = polynomials_over(mod7).factory

    # Test 1 root
    roots = [3]
    poly = zpoly(modulus, roots)
    assert poly(mod7(3)) == 0
    # Check equals x - 3
    assert poly == polysMod7([-3, 1])

    # Test 2 roots
    roots = [1, 2]
    poly = zpoly(modulus, roots)
    assert poly(mod7(1)) == 0
    assert poly(mod7(2)) == 0
    # Check equals x^2 - 3x + 2
    assert poly == polysMod7([2, -3, 1])


  def test_multi_inv(self):
    """Test of faster multiple inverse method."""
    # 6^-1 = 6
    modulus = 7
    mod7 = IntegersModP(modulus)
    outs = multi_inv(mod7, [mod7(6), mod7(6), mod7(6)])
    assert outs == [6, 6, 6]

    # 1^-1 = 1
    outs = multi_inv(mod7, [mod7(6), mod7(1), mod7(6)])
    assert outs == [6, 1, 6]

    outs = multi_inv(mod7, [mod7(0), mod7(1), mod7(1)])

    modulus = 2**256 - 2**32 * 351 + 1
    mod = IntegersModP(modulus)
    ## Root of unity such that x^precision=1
    G2 = mod(7)**((modulus - 1) // 4096)
    ### Powers of the higher-order root of unity
    xs = get_power_cycle(G2, modulus)
    xs_minus_1 = [x - 1 for x in xs]
    xs_minus_1_inv = multi_inv(mod, xs_minus_1)
    # Skip 0 since xs_minus_1[0] == 0
    for i in range(1, 5):
      assert xs_minus_1[i] * xs_minus_1_inv[i] == 1

    steps = 512
    precision = 4096
    z_evals = [xs[(i * steps) % precision] - 1 for i in range(precision)]
    z_inv = multi_inv(mod, z_evals)
    for i in range(1, 5):
      assert z_evals[i] * z_inv[i] == 1

  def test_lagrange_interp(self):
    """Test lagrangian interpolation."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysOverMod = polynomials_over(mod7).factory
    xs = [mod7(1), mod7(6)]
    ys = [mod7(1), mod7(6)]
    interp = lagrange_interp(modulus, xs, ys)
    # interp should equal x
    assert interp == polysOverMod([0, 1])

    xs = [mod7(1), mod7(6)]
    ys = [mod7(0), mod7(0)]
    interp = lagrange_interp(modulus, xs, ys)
    # interp should equal 0
    assert interp == polysOverMod([0])

  def test_lagrange_interp_4(self):
    """Test fast interpolation for degree 4 polynomials."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysOverMod = polynomials_over(mod7).factory
    xs = [mod7(1), mod7(2), mod7(3), mod7(6)]
    ys = [mod7(1), mod7(2), mod7(3), mod7(6)]
    interp = lagrange_interp_4(modulus, xs, ys)
    # interp should equal x
    assert interp == polysOverMod([0, 1])

  def test_lagrange_interp_2(self):
    """Test fast interpolation for degree 2 polynomials."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysOverMod = polynomials_over(mod7).factory
    xs = [mod7(1), mod7(2)]
    ys = [mod7(1), mod7(2)]
    interp = lagrange_interp_2(modulus, xs, ys)
    # interp should equal x
    assert interp == polysOverMod([0, 1])

  def test_multi_interp_4(self):
    """Test fast interpolation for multiple degree 4 polynomials."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysOverMod = polynomials_over(mod7).factory
    xs = [mod7(1), mod7(2), mod7(3), mod7(6)]
    ys = [mod7(1), mod7(2), mod7(3), mod7(6)]
    interp = multi_interp_4(modulus, [xs, xs], [ys, ys])
    # interp should equal x
    assert len(interp) == 2
    assert interp[0] == polysOverMod([0, 1])
    assert interp[1] == polysOverMod([0, 1])

  #def test_add_polys(self):
  #  """Test that addition of polynomials works."""
  #  field7 = PrimeField(7)
  #  added = field7.add_polys([1, 2], [1, 2])
  #  assert added == [2, 4]
