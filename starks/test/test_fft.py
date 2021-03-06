import unittest
from starks.fft import Additive_FFT
from starks.fft import NonBinaryFFT
from starks.fft import mul_polys
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.finitefield import FiniteField

class TestFFT(unittest.TestCase):
  """
  Basic tests for fft implementation. 
  """

  def test_Taylor_Expansion(self):
    p = 2
    m = 4
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    coefficients = [Zp(0)] * 5
    coefficients[0] = Zp(1)
    coefficients[1] = Zp(1)
    coefficients[4] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    # 1 + x + x^3
    f = field(polysOver([1, 1, 0, 1]))
    # 1 + x^2 + x^3
    fp = field(polysOver([1, 0, 1, 1]))
    obj = Additive_FFT(field)
    V1, V2 = obj.Taylor_Expansion(f, f.poly.degree())
    V1p, V2p = obj.Taylor_Expansion(fp, fp.poly.degree())

    assert V1.poly.degree() <= 3 and V2.poly.degree() <= 3
    assert V1p.poly.degree() <= 3 and V2p.poly.degree() <= 3

  def test_adfft(self):
    p = 2
    m = 4
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    coefficients = [Zp(0)] * 5
    coefficients[0] = Zp(1)
    coefficients[1] = Zp(1)
    coefficients[4] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    # 1 + x + x^3
    f = field(polysOver([1, 1, 0, 1]))
    obj = Additive_FFT(field)
    mp = 2
    beta = []
    beta.append(field(polysOver([1, 0, 1])))
    beta.append(field(polysOver([0, 1, 1])))
    shift = field(polysOver([1, 0, 1]))
    obj = Additive_FFT(field)
    V1, V2 = obj.Taylor_Expansion(f, f.poly.degree())
    W = obj.adfft(f, mp, beta)
    assert V1.poly.degree() <= 3 and V2.poly.degree() <= 3
    print(W)

  
  def test_adfft_inverse(self):
    p = 2
    m = 4
    Zp = IntegersModP(p)
    polysOver = polynomials_over(Zp)
    coefficients = [Zp(0)] * 5
    coefficients[0] = Zp(1)
    coefficients[1] = Zp(1)
    coefficients[4] = Zp(1)
    poly = polysOver(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)

    mp = 3
    x = []
    y = []
    x.append(field(polysOver([1, 0, 0])))
    y.append(field(polysOver([1, 1])))
    x.append(field(polysOver([1, 1, 1])))
    y.append(field(polysOver([0, 0, 1])))
    x.append(field(polysOver([1, 0])))
    y.append(field(polysOver([0, 1])))
    x.append(field(polysOver([1, 1])))
    y.append(field(polysOver([0, 1, 1])))
    x.append(field(polysOver([1, 1, 1])))
    y.append(field(polysOver([1, 0, 1])))
    x.append(field(polysOver([1, 0, 0])))
    y.append(field(polysOver([1, 1, 1])))
    x.append(field(polysOver([1])))
    y.append(field(polysOver([0])))
    x.append(field(polysOver([1])))
    y.append(field(polysOver([1, 1])))
    
    obj = Additive_FFT(field)
    f = obj.adfft_inverse(x, y, mp)
  

  def test_basic(self):
    """Basic test of fft."""
    modulus = 31 
    field = IntegersModP(31)
    polysOver = polynomials_over(field).factory
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = polysOver([val for val in range(4)])
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = field(3)**((modulus-1)//6)
    fft_solver = NonBinaryFFT(field, root_of_unity)
    evaluations = fft_solver.fft(poly)
    assert len(evaluations) == 6

  def test_large_modulus(self):
    """Basic test of fft with large modulus."""
    modulus = 2**256 - 2**32 * 351 + 1
    field = IntegersModP(modulus)
    polysOver = polynomials_over(field).factory
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = polysOver([val for val in range(4)])
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = field(7)**((modulus-1)//8)
    fft_solver = NonBinaryFFT(field, root_of_unity)
    evaluations = fft_solver.fft(poly)
    assert len(evaluations) == 8

  def test_fft_inv(self):
    """Test of Inverse FFT."""
    modulus = 31 
    field = IntegersModP(31)
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    polysOver = polynomials_over(field).factory
    poly = polysOver([val for val in range(4)])
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = field(3)**((modulus-1)//6)
    fft_solver = NonBinaryFFT(field, root_of_unity)
    evaluations = fft_solver.fft(poly)
    inv = fft_solver.inv_fft(evaluations)
    # Check we recover the original polynomial 
    assert inv == poly 

  def test_fft_output_type(self):
    """The output of FFT should be in the field if input is in field."""
    modulus = 31 
    field = IntegersModP(31)
    polysOver = polynomials_over(field).factory
    # 1 + 2x + 3x^2 + 4 x^3 mod 31
    poly = polysOver([val for val in range(4)])
    # TODO(rbharath): How does the choice of the n-th root of
    # unity make a difference in the fft?

    # A root of unity is a number such that z^n = 1
    # This provides us a 6-th root of unity (z^6 = 1)
    root_of_unity = field(3)**((modulus-1)//6)
    fft_solver = NonBinaryFFT(field, root_of_unity)
    evaluations = fft_solver.fft(poly)
    assert len(evaluations) == 6
    for val in evaluations:
      assert isinstance(val, field)

  # TODO(rbharath): Remove in future PR if confirmed not necessary
  #def test_fft_multidim(self):
  #  """Test FFT of multidimensional signal."""
  #  steps = 512 
  #  modulus = 2**256 - 2**32 * 351 + 1
  #  mod = IntegersModP(modulus)
  #  # [1 + 2x + 3x^2 + 4 x^3 mod 31,
  #  #  1 + 2x + 3x^2 + 4 x^3 mod 31]
  #  poly = [[mod(i), mod(i)] for i in range(4)]
  #  # Root of unity such that x^512=1
  #  G = mod(7)**((modulus - 1) // steps)
  #  evaluations = fft(poly, modulus, G, dims=2)
  #  assert len(evaluations) == steps
  #  assert len(evaluations[0]) == 2

  def test_mul_polys(self):
    """Test multiplication of polynomials."""
    modulus = 2**256 - 2**32 * 351 + 1
    mod = IntegersModP(modulus)
    # Root of unity such that x^512=1
    root_of_unity = mod(7)**((modulus - 1) // 512)
    # a = b = 1 + 2x + 3x^2 + 4 x^3 mod 31
    a = [mod(val) for val in range(4)] 
    b = [mod(val) for val in range(4)] 
    prod = mul_polys(a, b, root_of_unity)
