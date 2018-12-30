"""An implementation of Q/p

Q/p is the field of fractiosn over Z/p. This is canonically isomorphic to Z/p
but can be useful to have a separate representation for experimentation.
"""
from starks.euclidean import gcd
from starks.numbertype import FieldElement
from starks.numbertype import memoize
from starks.numbertype import typecheck
from starks.modp import _Modular


@memoize
def RationalsModP(p):
  """Assume p is a prime."""

  class RationalModP(_Modular):
    """A rational modulo p
    
    The rational is stored with numerator and denominator relatively prime.
    This is done to prevent growth in numerator or denominator which causes
    overflow and makes math harder. 
    """

    def __init__(self, m, n=1):
      """Constructor for rationals.

      If denominator is not specified, set it to 1.
      """
      try:
        # This is awkward constructor overloading
        if isinstance(m, bytes):
          num = int.from_bytes(m[:32], 'big')
          den = int.from_bytes(m[32:], 'big')
        else:
          num = int(m) % RationalModP.p
          den = int(n) % RationalModP.p
        common = gcd(num, den)
        # Handle case with 0
        if common == 0:
          self.m, self.n = num, den
        else:
          self.m = num // common
          self.n = den // common
      except:
        raise TypeError("Can't cast type %s to %s in __init__" %
                        (type(n).__name__, type(self).__name__))

      self.field = RationalModP

    @typecheck
    def __add__(self, other):
      num = (self.m * other.n + other.m * self.n) % RationalModP.p
      den = (self.n * other.n) % RationalModP.p
      common = gcd(num, den)
      return RationalModP(num // common, den // common)

    @typecheck
    def __sub__(self, other):
      num = (self.m * other.n - other.m * self.n) % RationalModP.p
      den = (self.n * other.n) % RationalModP.p
      common = gcd(num, den)
      return RationalModP(num // common, den // common)

    @typecheck
    def __mul__(self, other):
      num = (self.m * other.m) % RationalModP.p
      den = (self.n * other.n) % RationalModP.p
      return RationalModP(num, den)

    def __neg__(self):
      return RationalModP(-self.m, self.n)

    @typecheck
    def __eq__(self, other):
      return isinstance(other, RationalModP) and (
          (self.m * other.n) % RationalModP.p == (other.m * self.n) % RationalModP.p)

    @typecheck
    def __ne__(self, other):
      return isinstance(other, IntegerModP) is False or (
          (self.m * other.n) % RationalModP.p != (other.m * self.n) % RationalModP.p)

    # TODO(rbharath): This should be possible to implement. Think more about it.
    #@typecheck
    #def __divmod__(self, divisor):
    #  q, r = divmod(self.n, divisor.n)
    #  return (IntegerModP(q), IntegerModP(r))

    # TODO(rbharath): Check if this makes sense
    def inverse(self):
      if self.m == 0:
        raise Exception("Cannot invert with numerator 0")
      return RationalModP(self.n, self.m)
      ## need to use the division algorithm *as integers* because we're
      ## doing it on the modulus itself (which would otherwise be zero)
      #x, y, d = extended_euclidean_algorithm(self.n, self.p)

      #if d != 1:
      #  raise Exception("Error: p is not prime in %s!" % (self.__name__))

      #return IntegerModP(x)

    #def __abs__(self):
    #  return abs(self.n)

    def __str__(self):
      return "%s/%s" % (str(self.m), str(self.n))

    def __repr__(self):
      return '%d/%d (mod %d)' % (self.m, self.n, self.p)

    # TODO(rbharath): Can this method be done better?
    def to_bytes(self):
      return self.m.to_bytes(32, 'big') + self.n.to_bytes(32, 'big')

    #def __int__(self):
    #  return self.n

  RationalModP.p = p
  RationalModP.__name__ = 'Q/%d' % (p)
  RationalModP.englishName = 'RationalsMod%d' % (p)
  return RationalModP
