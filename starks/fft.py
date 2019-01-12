from typing import List
from starks.numbertype import FieldElement
from starks.numbertype import Vector

# TODO(rbharath): The type signatures here don't account for multidimensional inputs! Should this be List[Vector] instead?

class FFT(object):
  """Abstract class that specifies a FFT solver."""

  def __init__(self, field):
    self.field = field
    # How is this computed? Should this be a constructor argument?
    self.root_of_unity = None

  def fft(self, poly: Poly) -> List[FieldElement]:
    """The FFT efficiently evaluates a polynomial on many field elements."""
    raise NotImplementedError

  def inv_fft(self, List[FieldElement]) -> Poly:
    """Converts a polynomial represented as evaluations on m points to coefficients."""
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

def fft(vals: List[Vector], modulus: int, root_of_unity: FieldElement,
    inv: bool =False, dims:int =1) -> List[Vector]:
  """Computes FFT for potentially multidimensional sequences"""
  fft_vals = []
  for dim in range(dims):
    vals_dim = [val[dim] for val in vals]
    fft_dim = fft_1d(vals_dim, modulus, root_of_unity, inv=inv)
    fft_vals.append(fft_dim)
  # We get tuples without the explicit list cast
  fft_joint = list([list(elt) for elt in zip(*fft_vals)])
  return fft_joint


def fft_1d(vals: List[FieldElement], modulus: int, root_of_unity: FieldElement, inv: bool = False) -> List[FieldElement]:
  """Computes FFT for one dimensional inputs"""
  # Build up roots of unity
  rootz = [1, root_of_unity]
  while rootz[-1] != 1:
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
