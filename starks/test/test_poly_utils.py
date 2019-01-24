import unittest
from starks.modp import IntegersModP
from starks.multivariate_polynomial import multivariates_over
from starks.poly_utils import zpoly
from starks.poly_utils import multi_inv
from starks.poly_utils import lagrange_interp
from starks.poly_utils import lagrange_interp_2
from starks.poly_utils import lagrange_interp_4
from starks.poly_utils import multi_interp_4 
from starks.poly_utils import is_primitive
from starks.poly_utils import is_monic
from starks.poly_utils import is_irreducible
from starks.poly_utils import generate_primitive_polynomial 
from starks.polynomial import polynomials_over
from starks.poly_utils import construct_multivariate_dirac_delta
from starks.poly_utils import construct_multivariate_coefficients
from starks.poly_utils import project_to_univariate
from starks.utils import get_power_cycle
from starks.finitefield import FiniteField

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

  def test_project_multivar(self):
    """Test the multivariate projection."""
    modulus = 7
    field = IntegersModP(modulus)
    width = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(field, width).factory
    # This should equal xy
    xy_poly = multi({(1, 1, 0): 1})
    x_poly = project_to_univariate(xy_poly, 0, field, width)
    # TODO(rbharath): Add more nontrivial tests

  def test_zpoly(self):
    """Test construction of polynomials with specified root"""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysMod7 = polynomials_over(mod7).factory

    # Test 1 root
    roots = [3]
    poly = zpoly(mod7, roots)
    assert poly(mod7(3)) == 0
    # Check equals x - 3
    assert poly == polysMod7([-3, 1])

    # Test 2 roots
    roots = [1, 2]
    poly = zpoly(mod7, roots)
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
    field = IntegersModP(modulus)
    ## Root of unity such that x^precision=1
    G2 = field(7)**((modulus - 1) // 4096)
    ### Powers of the higher-order root of unity
    xs = get_power_cycle(G2, field)
    xs_minus_1 = [x - 1 for x in xs]
    xs_minus_1_inv = multi_inv(field, xs_minus_1)
    # Skip 0 since xs_minus_1[0] == 0
    for i in range(1, 5):
      assert xs_minus_1[i] * xs_minus_1_inv[i] == 1

    steps = 512
    precision = 4096
    z_evals = [xs[(i * steps) % precision] - 1 for i in range(precision)]
    z_inv = multi_inv(field, z_evals)
    for i in range(1, 5):
      assert z_evals[i] * z_inv[i] == 1

  def test_lagrange_interp(self):
    """Test lagrangian interpolation."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysOverMod = polynomials_over(mod7).factory
    xs = [mod7(1), mod7(6)]
    ys = [mod7(1), mod7(6)]
    interp = lagrange_interp(mod7, xs, ys)
    # interp should equal x
    assert interp == polysOverMod([0, 1])

    xs = [mod7(1), mod7(6)]
    ys = [mod7(0), mod7(0)]
    interp = lagrange_interp(mod7, xs, ys)
    # interp should equal 0
    assert interp == polysOverMod([0])

    # Test lagrange interp over general finite field
    Z5 = IntegersModP(5)
    F25 = FiniteField(5, 2)
    polysOverF = polynomials_over(F25).factory
    xs = [F25(1), F25(2)]
    ys = [F25(1), F25(2)]
    interp = lagrange_interp(F25, xs, ys)
    # interp should equal x
    assert interp == polysOverF([F25(0), F25(1)])

  def test_lagrange_interp_4(self):
    """Test fast interpolation for degree 4 polynomials."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysOverMod = polynomials_over(mod7).factory
    xs = [mod7(1), mod7(2), mod7(3), mod7(6)]
    ys = [mod7(1), mod7(2), mod7(3), mod7(6)]
    #interp = lagrange_interp_4(modulus, xs, ys)
    interp = lagrange_interp_4(mod7, xs, ys)
    # interp should equal x
    assert interp == polysOverMod([0, 1])

  def test_lagrange_interp_2(self):
    """Test fast interpolation for degree 2 polynomials."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysOverMod = polynomials_over(mod7).factory
    xs = [mod7(1), mod7(2)]
    ys = [mod7(1), mod7(2)]
    #interp = lagrange_interp_2(modulus, xs, ys)
    interp = lagrange_interp_2(mod7, xs, ys)
    # interp should equal x
    assert interp == polysOverMod([0, 1])

  def test_multi_interp_4(self):
    """Test fast interpolation for multiple degree 4 polynomials."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    polysOverMod = polynomials_over(mod7).factory
    xs = [mod7(1), mod7(2), mod7(3), mod7(6)]
    ys = [mod7(1), mod7(2), mod7(3), mod7(6)]
    #interp = multi_interp_4(modulus, [xs, xs], [ys, ys])
    interp = multi_interp_4(mod7, [xs, xs], [ys, ys])
    # interp should equal x
    assert len(interp) == 2
    assert interp[0] == polysOverMod([0, 1])
    assert interp[1] == polysOverMod([0, 1])

  def test_is_primitive(self):
    """Tests whether the primitivity check is correctly implemented."""
    modulus = 2
    degree = 2
    mod = IntegersModP(modulus)
    polysOver = polynomials_over(mod).factory

    # From table 4.6 in
    # http://math.fau.edu/bkhadka/Syllabi/A%20handbook%20of%20applied%20cryptography.pdf
    # x^2 + x + 1 is primitive over Z/2
    prim_poly = polysOver([1, 1, 1])
    assert is_primitive(prim_poly, modulus, degree)

    # x^2 is not primitive over Z/2 
    x_square = polysOver([0, 0, 1])
    assert not is_primitive(x_square, modulus, degree)

    # x^9 + x + 1 is primitive over Z/2
    coeffs = [0] * 10
    coeffs[0] = 1
    coeffs[1] = 1
    coeffs[-1] = 1
    prim_poly = polysOver(coeffs)
    assert is_primitive(prim_poly, modulus, degree)

  def test_generate_primitive_poly(self):
    """Tests the generation of primitive polynomials."""
    modulus = 2
    degree = 2
    gen_poly = generate_primitive_polynomial(modulus, degree)
    assert is_irreducible(gen_poly, modulus)
    assert is_primitive(gen_poly, modulus, degree)

    degree = 5
    gen_poly = generate_primitive_polynomial(modulus, degree)
    assert is_irreducible(gen_poly, modulus)
    assert is_primitive(gen_poly, modulus, degree)

  def test_is_monic(self):
    """Tests the is_monic primitive."""
    modulus = 3
    degree = 2
    mod = IntegersModP(modulus)
    polysOver = polynomials_over(mod).factory
    # x^2 + x + 1 is monic over Z/2
    poly = polysOver([1, 1, 1])
    assert is_monic(poly)

    # 2x^2 + x + 1 is not monic over Z/2
    poly = polysOver([1, 1, 2])
    assert not is_monic(poly)

  def test_construct_multivariate_dirac_delta(self):
    """Tests the construct of the multivariate dirac delta."""
    modulus = 3
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory
    # Let's generate the dirac delta at x=0, y=0, z=0
    values = [mod7(0), mod7(0), mod7(0)]
    dirac = construct_multivariate_dirac_delta(mod7, values, n)

    # The dirac delta should be 1 at x=0, y=0, z=0
    assert dirac((0, 0, 0)) == 1
    # It should be 0 elsewhere
    assert dirac((1, 0, 0)) == 0
    assert dirac((0, 1, 0)) == 0
    assert dirac((0, 0, 1)) == 0

    # Let's generate the dirac delta at x=1, y=1, z=1
    values = [mod7(1), mod7(1), mod7(1)]
    dirac = construct_multivariate_dirac_delta(mod7, values, n)

    # The dirac delta should be 1 at x=1, y=1, z=1
    assert dirac((1, 1, 1)) == 1
    # It should be 0 elsewehre
    assert dirac((1, 0, 0)) == 0
    assert dirac((0, 1, 0)) == 0
    assert dirac((0, 0, 1)) == 0

  def test_construct_multivariate_coefficients(self):
    """Tests the "compilation" of a function into a polynomial."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    n = 3
    # Let's make polynomials in (Z/7)[x, y, z]
    multi = multivariates_over(mod7, n).factory
    # Our test function
    def f(x):
      if isinstance(x, tuple) or isinstance(x, list):
        x = x[0]
      return x + 1
    poly = construct_multivariate_coefficients(mod7, f, n)

    # This should equal x + 1 
    x_plus_one_poly = multi({(0, 0, 0): mod7(1), (1, 0, 0): mod7(1)})
    assert poly == x_plus_one_poly

