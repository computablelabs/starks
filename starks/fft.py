from typing import List
from starks.numbertype import Field
from starks.numbertype import FieldElement
from starks.numbertype import Vector
from starks.numbertype import Poly
from starks.polynomial import polynomials_over
from starks.modp import IntegersModP

def int_to_bin_string(i):
    if i == 0:
        return "0"
    s = ''
    while i:
        if i & 1 == 1:
            s = "1" + s
        else:
            s = "0" + s
        i //= 2
    return s

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


class Additive_FFT(FFT):
  def __init__(self, field):
    self.field = field

  def Taylor_Expansion(self, Polys, n):
    if n <= 2:
      return Polys

    for x in range(n):
      if 2**(x+1) < n and 2**(x+2) >= n:
        k = x

    polysOver = polynomials_over(IntegersModP(2))
    list_f0 = []
    for i in range(2**(k+1)):
      if i > Polys.poly.degree():
        list_f0.append(0)
      else:
        list_f0.append(int(str(Polys.poly.coefficients[i])[0]))

    list_f1 = []
    for i in range(2**(k)):
      if 2**(k+1)+i > Polys.poly.degree():
        list_f1.append(0)
      else:
        list_f1.append(int(str(Polys.poly.coefficients[2**(k+1)+i])[0]))

    list_f2 = []
    for i in range(2**k):
      if 2**(k+1)+2**k+i > Polys.poly.degree():
        list_f2.append(0)
      else:
        list_f2.append(int(str(Polys.poly.coefficients[2**(k+1)+2**k+i])[0]))

    f0 = self.field(polysOver(list_f0))
    f1 = self.field(polysOver(list_f1))
    f2 = self.field(polysOver(list_f2))

    h = f1+f2
    twoK = []
    for i in range(2**(k)):
      twoK.append(0)
    twoK.append(1)
    f_twoK = self.field(polysOver(twoK))

    g0 = f0+f_twoK*h
    g1 = h+f_twoK*f2

    V1 = self.Taylor_Expansion(g0, n/2)
    V2 = self.Taylor_Expansion(g1, n/2)

    return V1, V2



  def adfft(self, Polys, m, affine_beta):
    polysOver = polynomials_over(IntegersModP(2))
    f1 = self.field(polysOver([0]))
    for i in range(Polys.poly.degree()+1):
      if str(Polys.poly.coefficients[i])[0] == '1':
        f1 = f1 + (affine_beta[0])**i

    f2 = self.field(polysOver([0]))
    for i in range(Polys.poly.degree()+1):
      if str(Polys.poly.coefficients[i])[0] == '1':
        f2 = f2

    if m == 1:
      return f1, f2

    g = self.field(polysOver([0]))
    x = self.field(polysOver([0, 1]))
    x = x * affine_beta[m-1] 
    for i in range(Polys.poly.degree()+1):
      if str(Polys.poly.coefficients[i])[0] == '1':
        g = g + x**i

    g0, g1 = self.Taylor_Expansion(g, 2**m)

    gamma = []
    beta_m_I = affine_beta[m-1].inverse()
    for i in range(m-1):
      gamma.append(affine_beta[i] * beta_m_I);

    delta = []
    for i in range(m-1):
      delta.append(gamma[i]**2-gamma[i])

    G = []
    for i in range(2**(m-1)):
      binary = int_to_bin_string(i)
      temp = self.field(polysOver([0]))
      for j in range(len(binary)):
        if binary[j] == '1':
          temp = temp + gamma[j]
      G.append(temp)

    D = delta

    u = self.adfft(g0, m-1, D)
    v = self.adfft(g1, m-1, D)

    w1 = []
    w2 = []
    for i in range(2**(m-1)):
      w1.append(u[i]+G[i]*v[i])
      w2.append(w1[i]+v[i])

    w = []
    for i in range(len(w1)):
      w.append(w1[i])
    for i in range(len(w2)):
      w.append(w2[i])

    return w

  def adfft_inverse(self, x, y, m): 
    if m == 1:
      if x[0] == x[1]:
        return x[0]
      else:
        return ((y[1]-y[0])/(x[1]-x[0]))*(self.field(polysOver([0, 1]))-x[0])+y[0]

    polysOver = polynomials_over(IntegersModP(2))
    beta = []
    for i in range(m):
      beta.append(x[2**i])

    gamma = []
    Ibeta = beta[-1].inverse()
    for i in range(m-1):
      gamma.append(beta[i] * Ibeta);

    delta = []
    for i in range(m-1):
      delta.append(gamma[i]*gamma[i]-gamma[i])

    G = []
    for i in range(2**(m-1)):
      binary = int_to_bin_string(i)
      temp = self.field(polysOver([0]))
      for j in range(len(binary)):
        if binary[j] == '1':
          temp = temp + gamma[j]
      G.append(temp)

    D = delta

    v = []
    u = []
    for i in range(2**(m-1)):
      v.append(y[i+2**(m-1)]-y[i])
      u.append(y[i] - G[i]*v[i])

    x1 = []
    for i in range(2**(m-1)):
      binary = int_to_bin_string(i+1)
      temp = self.field(polysOver([0]))
      for j in range(len(binary)):
        if binary[j] == '1':
          temp = temp + D[j]
      x1.append(temp)

    g_0 = self.adfft_inverse(x1, u, m-1)
    g_1 = self.adfft_inverse(x1, v, m-1)

    g = self.field(polysOver([0]))
    g_right_tempp = self.field(polysOver([0, 1])) * Ibeta
    g_right_temp = g_right_tempp * g_right_tempp - g_right_tempp
    g_right = []
    g_right.append(self.field(polysOver([1])))
    multiplier = []
    multiplier.append(self.field(polysOver([0]))) 
    multiplier.append(g_right_tempp)
    multiplier.append(self.field(polysOver([1])))
    multiplier.append(self.field(polysOver([1]))+g_right_tempp)
    for i in range(2**(m-1)):
      g_right.append(g_right[-1]*g_right_temp)

    for i in range(2**(m-1)):
      if i <= g_0.poly.degree():
        g0I = int(g_0.poly.coefficients[i])
      else:
        g0I = 0

      if i <= g_1.poly.degree():
        g1I = int(g_1.poly.coefficients[i])
      else:
        g1I = 0
      g  = g + multiplier[g0I+2**g1I] * g_right[i]

    return g


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
