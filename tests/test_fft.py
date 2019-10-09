import pytest
#from starks.fft import NonBinaryFFT
#from starks.fft import mul_polys
#from starks.modp import IntegersModP
#from starks.polynomial import polynomials_over

#class TestFFT(unittest.TestCase):
#  """
#  Basic tests for fft implementation. 
#  """
#  def test_Taylor_Expansion(self):
#    polysOver = polynomials_over(IntegersModP(2))
#    # 1 + x + x^3
#    f = field(polysOver([1, 1, 0, 1]))
#    V1, V2 = Taylor_Expansion(f, f.degree())
#
#    print(f)
#    print(V1)
#    print(V2)
#
#
#  def test_adfft(self):
#    polysOver = polynomials_over(IntegersModP(2))
#    # 1 + x + x^3 + x^7
#    f = field(polysOver([1, 1, 0, 1, 0, 0, 0, 1]))
#    m = 3
#    beta = []
#    beta.append(field(polysOver([1, 0, 0, 1, 1])))
#    beta.append(field(polysOver([1, 1, 0, 1])))
#    beta.append(field(polysOver([1, 0, 1, 1])))
#    shift = field(polysOver([1, 0, 0, 0, 1]))
#
#    V1, V2 = Taylor_Expansion(f, f.degree())
#    print(V1)
#    print(V2)
#
#    W = adfft(f, m, beta, shift)
#    print(W)
#  
#  def test_basic(self):
#    """Basic test of fft."""
#    modulus = 31 
#    field = IntegersModP(31)
#    polysOver = polynomials_over(field).factory
#    # 1 + 2x + 3x^2 + 4 x^3 mod 31
#    poly = polysOver([val for val in range(4)])
#    # TODO(rbharath): How does the choice of the n-th root of
#    # unity make a difference in the fft?
#
#    # A root of unity is a number such that z^n = 1
#    # This provides us a 6-th root of unity (z^6 = 1)
#    root_of_unity = field(3)**((modulus-1)//6)
#    fft_solver = NonBinaryFFT(field, root_of_unity)
#    evaluations = fft_solver.fft(poly)
#    assert len(evaluations) == 6
#
#  def test_large_modulus(self):
#    """Basic test of fft with large modulus."""
#    modulus = 2**256 - 2**32 * 351 + 1
#    field = IntegersModP(modulus)
#    polysOver = polynomials_over(field).factory
#    # 1 + 2x + 3x^2 + 4 x^3 mod 31
#    poly = polysOver([val for val in range(4)])
#    # TODO(rbharath): How does the choice of the n-th root of
#    # unity make a difference in the fft?
#
#    # A root of unity is a number such that z^n = 1
#    # This provides us a 6-th root of unity (z^6 = 1)
#    root_of_unity = field(7)**((modulus-1)//8)
#    fft_solver = NonBinaryFFT(field, root_of_unity)
#    evaluations = fft_solver.fft(poly)
#    assert len(evaluations) == 8
#
#  def test_fft_inv(self):
#    """Test of Inverse FFT."""
#    modulus = 31 
#    field = IntegersModP(31)
#    # 1 + 2x + 3x^2 + 4 x^3 mod 31
#    polysOver = polynomials_over(field).factory
#    poly = polysOver([val for val in range(4)])
#    # TODO(rbharath): How does the choice of the n-th root of
#    # unity make a difference in the fft?
#
#    # A root of unity is a number such that z^n = 1
#    # This provides us a 6-th root of unity (z^6 = 1)
#    root_of_unity = field(3)**((modulus-1)//6)
#    fft_solver = NonBinaryFFT(field, root_of_unity)
#    evaluations = fft_solver.fft(poly)
#    inv = fft_solver.inv_fft(evaluations)
#    # Check we recover the original polynomial 
#    assert inv == poly 
#
#  def test_fft_output_type(self):
#    """The output of FFT should be in the field if input is in field."""
#    modulus = 31 
#    field = IntegersModP(31)
#    polysOver = polynomials_over(field).factory
#    # 1 + 2x + 3x^2 + 4 x^3 mod 31
#    poly = polysOver([val for val in range(4)])
#    # TODO(rbharath): How does the choice of the n-th root of
#    # unity make a difference in the fft?
#
#    # A root of unity is a number such that z^n = 1
#    # This provides us a 6-th root of unity (z^6 = 1)
#    root_of_unity = field(3)**((modulus-1)//6)
#    fft_solver = NonBinaryFFT(field, root_of_unity)
#    evaluations = fft_solver.fft(poly)
#    assert len(evaluations) == 6
#    for val in evaluations:
#      assert isinstance(val, field)
#
#  # TODO(rbharath): Remove in future PR if confirmed not necessary
#  #def test_fft_multidim(self):
#  #  """Test FFT of multidimensional signal."""
#  #  steps = 512 
#  #  modulus = 2**256 - 2**32 * 351 + 1
#  #  mod = IntegersModP(modulus)
#  #  # [1 + 2x + 3x^2 + 4 x^3 mod 31,
#  #  #  1 + 2x + 3x^2 + 4 x^3 mod 31]
#  #  poly = [[mod(i), mod(i)] for i in range(4)]
#  #  # Root of unity such that x^512=1
#  #  G = mod(7)**((modulus - 1) // steps)
#  #  evaluations = fft(poly, modulus, G, dims=2)
#  #  assert len(evaluations) == steps
#  #  assert len(evaluations[0]) == 2
#
#  def test_mul_polys(self):
#    """Test multiplication of polynomials."""
#    modulus = 2**256 - 2**32 * 351 + 1
#    mod = IntegersModP(modulus)
#    # Root of unity such that x^512=1
#    root_of_unity = mod(7)**((modulus - 1) // 512)
#    # a = b = 1 + 2x + 3x^2 + 4 x^3 mod 31
#    a = [mod(val) for val in range(4)] 
#    b = [mod(val) for val in range(4)] 
#    prod = mul_polys(a, b, root_of_unity)
#
