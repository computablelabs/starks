from typing import List
from starks.numbertype import Field
from starks.numbertype import FieldElement
from starks.numbertype import Vector
from starks.numbertype import Poly
from starks.polynomial import polynomials_over

class Additive_FFT(object):
  def __init__(self):
    raise NotImplementedError

  def Taylor_Expansion(self, Polys, degree):
    if degree <= 2:
      return Polys

    for x in range(degree):
      if 2**(x+1) < degree and 2**(x+2) >= degree:
        k = x

    polysOver = polynomials_over(IntegersModP(2))
    list_f0 = []
    for i in range(2**(k+1)):
      list_f0.append(Polys.coefficients[i])

    list_f1 = []
    for i in range(2**(k)):
        list_f0.append(Polys.coefficients[2**(k+1)+i])

    list_f2 = []
    for i in range(2**k):
        list_f0.append(Polys.coefficients[2**(k+1)+2**k+i])

    f0 = field(polysOver(list_f0))
    f1 = field(polysOver(list_f1))
    f2 = field(polysOver(list_f2))

    h = f1+f2
    twoK = []
    for i in range(2**(k)):
      twoK.append(0)
      twoK.append(1)
      f_twoK = field(polysOver(twoK))

    g0 = f0+f_twoK*h
    g1 = h+f_twoK*f2

    V1 = Taylor_Expansion(g0, degree/2)
    V2 = Taylor_Expansion(g1, degree/2)

    return V1, V2



def adfft(self, Polys, m, affine_beta, shift):
  if m == 1:
    return Polys(shift), Polys(shift+affine_beta[0])

  polysOver = polynomials_over(IntegersModP(2))
  list_g = []
  for i in range(Polys.degree()+1):
    list_g.append((affine_beta[m-1]**i*Polys.coefficients[i])%2)

    g = field(polysOver(list_g))

    g0, g1 = Taylor_Expansion(g, g.degree())

    gamma = []
    for i in range(m-1):
      gamma.append(affine_beta[i] * affine_beta[m-1]).inverse();

    delta = []
    for i in range(m-1):
      delta.append(gamma[i]**2-gamma[i])

    S_G = shift * affine_beta[m-1].inverse()
    S_D = S_G**2-S_G

    G = []
    G.append(S_G)
    for i in range(m-1):
      G.append(gamma[i])

    D = delta

    u = Additive_FFT(g0, m-1, D, S_D)
    v = Additive_FFT(g1, m-1, D, S_D)

    w1 = []
    for i in range(2**(m-1)):
      w1.append(u[i]+G[i]*u[i])
      w2.append(w[i]+v[i])

    w = w1 + w2
    return w

def adfft_inverse(self, x, y, m): 
  beta = []
  for i in range(m):
    beta.append(x[2**i])

  gamma = []
  for i in range(m-1):
    gamma.append(beta[i] * beta[m-1]).inverse();

  delta = []
  for i in range(m-1):
    delta.append(gamma[i]**2-gamma[i])

  S_G = shift * affine_beta[m-1].inverse()
  S_D = S_G**2-S_G

  G = []
  G.append(S_G)
  for i in range(m-1):
    G.append(gamma[i])

  D = delta

  v = []
  u = []
  for i in range(2**(m-1)):
    v.append(y[i+2**(m-1)]-w[i])
    u.append(w[i] - G[i]*v[i])

  g_0 = self.adfft_inverse(D, u, m-1)
  g_1 = self.adfft_inverse(D, v, m-1)

  g = field(polysOver([0]))
  g_right_temp = field(polysOver([0, 1, 1]))
  g_right = []
  g_right.append(field(polysOver([1])))
  multiplier = []
  multiplier.append(field(polysOver([0]))) 
  multiplier.append(field(polysOver([0, 1])))
  multiplier.append(field(polysOver([1])))
  multiplier.append(field(polysOver([1, 1])))
  for i in range(2**(m-1)-1):
    g_right.append(g_right[-1]*g_right_temp)

  for i in range(2**(m-1)):
    g  = g + multiplier[g_0.coefficients[i]+2**g_1.coefficients[i]] * g_right[i]

  return g


class FFT(object):
  """Abstract class that specifies a FFT solver."""

  def __init__(self):
    raise NotImplementedError

  def fft(self, poly: Poly) -> List[FieldElement]:
    """The FFT efficiently evaluates a polynomial on many field elements."""
    raise NotImplementedError

  def inv_fft(self, values: List[FieldElement]) -> Poly:
    """Converts a polynomial represented as evaluations on m points to coefficients."""
    raise NotImplementedError


class NonBinaryFFT(FFT):
  """FFT that works for finite fields which don't have characteristic 2."""
  def __init__(self, field, root_of_unity):
    self.field = field
    self.root_of_unity = root_of_unity
    self.polysOver = polynomials_over(self.field).factory

  def fft(self, poly: Poly) -> List[FieldElement]:
    """Runs FFT algorithm."""
    return fft_1d(self.field, poly.coefficients, self.field.p, self.root_of_unity,
        inv=False)

  def inv_fft(self, values: List[FieldElement]) -> Poly:
    """Performs the inverse fft."""
    coeffs = fft_1d(self.field, values, self.field.p, self.root_of_unity,
        inv=True)
    return self.polysOver(coeffs)

class BinaryFFT(FFT):
  """FFT that works for finite fields of characteristic 2.

  Implements basis and algorithm from https://arxiv.org/pdf/1404.3458.pdf"""

  def fft(self, poly: Poly) -> List[FieldElement]:
    """Runs FFT algorithm."""
    raise NotImplementedError

  def inv_fft(self, values: List[FieldElement]) -> Poly:
    """Performs the inverse fft."""
    raise NotImplementedError

def _simple_ft(vals: List[FieldElement], roots_of_unity: FieldElement) -> List[FieldElement]:
  """Efficient base case implementation.
  
  The FFT recurses down halves of the list. This method is
  called to handle the base case of the fft.
  """
  L = len(roots_of_unity)
  o = []
  for i in range(L):
    last = 0
    for j in range(L):
      last += vals[j] * roots_of_unity[(i * j) % L]
    o.append(last)
  return o


def _fft(vals: List[FieldElement], roots_of_unity: FieldElement) -> List[FieldElement]:
  if len(vals) <= 4:
    #return vals
    return _simple_ft(vals, roots_of_unity)
  L = _fft(vals[::2], roots_of_unity[::2])
  R = _fft(vals[1::2], roots_of_unity[::2])
  o = [0 for i in vals]
  for i, (x, y) in enumerate(zip(L, R)):
    y_times_root = y * roots_of_unity[i]
    o[i] = (x + y_times_root)
    o[i + len(L)] = (x - y_times_root)
  return o

def fft_1d(field: Field, vals: List[FieldElement], modulus: int, root_of_unity: FieldElement, inv: bool = False) -> List[FieldElement]:
  """Computes FFT for one dimensional inputs"""
  # Build up roots of unity
  rootz = [field(1), root_of_unity]
  while rootz[-1] != field(1):
    rootz.append((rootz[-1] * root_of_unity))
  # Fill in vals with zeroes if needed
  if len(rootz) > len(vals) + 1:
    vals = vals + [0] * (len(rootz) - len(vals) - 1)
  if inv:
    # Inverse FFT
    invlen = pow(len(vals), modulus - 2, modulus)
    return [(x * invlen) for x in _fft(vals, rootz[:0:-1])]
  else:
    # Regular FFT
    return _fft(vals, rootz[:-1])


def mul_polys(a: List[FieldElement], b: List[FieldElement], root_of_unity: FieldElement) -> List[FieldElement]:
  """Multiply polynomials by converting to fourier space"""
  rootz = [1, root_of_unity]
  while rootz[-1] != 1:
    rootz.append((rootz[-1] * root_of_unity))
  if len(rootz) > len(a) + 1:
    a = a + [0] * (len(rootz) - len(a) - 1)
  if len(rootz) > len(b) + 1:
    b = b + [0] * (len(rootz) - len(b) - 1)
  x1 = _fft(a, rootz[:-1])
  x2 = _fft(b, rootz[:-1])
  return _fft([(v1 * v2) for v1, v2 in zip(x1, x2)], rootz[:0:-1])
