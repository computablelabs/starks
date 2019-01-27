"""An implementation of Z/p.

Adapted from https://github.com/j2kun/finite-fields
"""

from starks.euclidean import extended_euclidean_algorithm
from starks.numbertype import FieldElement
from starks.numbertype import memoize
from starks.numbertype import typecheck

#from starks.numbertype import *


# so all IntegersModP are instances of the same base class
class _Modular(FieldElement):
  pass

# TODO(rbharath): This sort of memoized architecture feels awkward. See if this
# is good style or not.
@memoize
def IntegersModP(p):
  """Assume p is prime"""

  class IntegerModP(_Modular):

    def __init__(self, n):
      try:
        if isinstance(n, bytes):
          self.n = int.from_bytes(n, 'big')
        else:
          ##############################################################
          try:
            a = int(n)
          except Exception:
            print("n")
            print(n)
          ##############################################################
          self.n = int(n) % IntegerModP.p
      except:
        raise TypeError("Can't cast type %s to %s in __init__" %
                        (type(n).__name__, type(self).__name__))

      self.field = IntegerModP

    @typecheck
    def __add__(self, other):
      return IntegerModP(self.n + other.n)

    @typecheck
    def __sub__(self, other):
      return IntegerModP(self.n - other.n)

    @typecheck
    def __mul__(self, other):
      return IntegerModP(self.n * other.n)

    def __neg__(self):
      return IntegerModP(-self.n)

    @typecheck
    def __eq__(self, other):
      return isinstance(other, IntegerModP) and self.n == other.n

    @typecheck
    def __ne__(self, other):
      return isinstance(other, IntegerModP) is False or self.n != other.n

    @typecheck
    def __divmod__(self, divisor):
      q, r = divmod(self.n, divisor.n)
      return (IntegerModP(q), IntegerModP(r))

    def inverse(self):
      # need to use the division algorithm *as integers* because we're
      # doing it on the modulus itself (which would otherwise be zero)
      x, y, d = extended_euclidean_algorithm(self.n, self.p)

      if d != 1:
        raise Exception("Error: p is not prime in %s!" % (self.__name__))

      return IntegerModP(x)

    def __abs__(self):
      return abs(self.n)

    def __str__(self):
      return str(self.n)

    def __repr__(self):
      return '%d (mod %d)' % (self.n, self.p)

    def __int__(self):
      return self.n

    # TODO(rbharath): Can this method be done better?
    def to_bytes(self):
      return self.n.to_bytes(32, 'big')


  IntegerModP.p = p
  # To match finite field parameters.
  IntegerModP.m = 1
  IntegerModP.__name__ = 'Z/%d' % (p)
  IntegerModP.englishName = 'IntegersMod%d' % (p)
  IntegerModP.field_size = p
  return IntegerModP
