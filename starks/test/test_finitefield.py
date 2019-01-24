import unittest
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.finitefield import FiniteField
from starks.finitefield import generate_irreducible_polynomial 
from starks.poly_utils import is_irreducible

class TestFiniteField(unittest.TestCase):
  """"
  Basic tests for finite field construction.
  """

  def test_construction(self):
    """Test constructors work."""
    F23 = FiniteField(2, 3)
    x = F23([1, 1])

    F35 = FiniteField(3, 5)
    y = F35([1, 1, 2])

  def test_generator(self):
    """Test the construction of a field generator."""
    p = 2
    m = 8
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    field = FiniteField(p, m)
    # X is a gen
    X = field(polysOver([0, 1]))
    order = 1
    val = X
    while val != 1:
      val *= X
      order += 1
    assert order == p**m - 1

  def test_large(self):
    """Tests construction of large finite field."""
    p = 2
    m = 17 
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    #field = FiniteField(p, m)
    #x^17 + x^3 + 1 is primitive 
    coefficients = [Zp(0)] * 18
    coefficients[0] = Zp(1)
    coefficients[3] = Zp(1)
    coefficients[17] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    # The fact this field can be constructed means an error
    # wasn't raised and that poly is primitive.


  def test_irreducibility(self):
    """Test the irreducibility algorithm"""
    def p(L, q):
      f = IntegersModP(q)
      Polynomial = polynomials_over(f).factory
      return Polynomial(L)
    # p(x) = x over Z/2
    assert is_irreducible(p([0,1], 2), 2)
    # p(x) = 1 + x^2 over Z/2
    assert not is_irreducible(p([1,0,1], 2), 2)
    # p(x) = 1 + x^2 over Z/3
    assert is_irreducible(p([1,0,1], 3), 3)

    # p(x) = 1 + x^3 over Z/5
    assert not is_irreducible(p([1,0,0,1], 5), 5)
    # p(x) = 1 + x^3 over Z/7
    assert not is_irreducible(p([1,0,0,1], 7), 7)
    # p(x) = 1 + x^3 over Z/11
    assert not is_irreducible(p([1,0,0,1], 11), 11)

    # p(x) = -2 + x^2 over Z/13
    assert is_irreducible(p([-2, 0, 1], 13), 13)

  def test_basic_finite_field(self):
    """Do some basic finite field tests."""
    Z5 = IntegersModP(5)
    Poly = polynomials_over(Z5).factory
    f = Poly([3,0,1])
    F25 = FiniteField(5, 2, polynomialModulus=f)
    x = F25([2,1])
    assert Poly([1,2]) == x.inverse()

  def test_generate_irreducible_polynomial(self):
    """Test generation of irreducible polynomials."""
    modulus = 2
    degree = 4 
    poly = generate_irreducible_polynomial(modulus, degree)
    assert is_irreducible(poly, modulus)

  def test_large_finite_field(self):
    """Test creation of GF(2^16)"""
    p = 2
    m = 16
    field = FiniteField(p, m)
