"""An implementation of the Euclidean algorithm

Adapted from https://github.com/j2kun/finite-fields
"""


def gcd(a, b):
  """General implementation of Euclidean algorithm.
  
  A general Euclidean algorithm for any number type with a divmod and a
  valuation abs() whose minimum value is zero. Note that this implementation
  can work with objects and not just native Python bigints.
  """
  if abs(a) < abs(b):
    return gcd(b, a)

  while abs(b) > 0:
    _, r = divmod(a, b)
    a, b = b, r

  return a


def extended_euclidean_algorithm(a, b):
  """extended_euclidean_algorithm: int, int -> int, int, int

  Input (a,b) and output three numbers x,y,d such that ax + by = d = gcd(a,b).
  Works for any number type with a divmod and a valuation abs()
  whose minimum value is zero
  """
  if abs(b) > abs(a):
    (x, y, d) = extended_euclidean_algorithm(b, a)
    return (y, x, d)

  if abs(b) == 0:
    return (1, 0, a)

  x1, x2, y1, y2 = 0, 1, 1, 0
  while abs(b) > 0:
    q, r = divmod(a, b)
    x = x2 - q * x1
    y = y2 - q * y1
    a, b, x2, x1, y2, y1 = b, r, x1, x, y1, y

  return (x2, y2, a)
