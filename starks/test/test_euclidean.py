import unittest
from starks.euclidean import gcd
from starks.euclidean import extended_euclidean_algorithm


def run_foo(expected, actual):
  if expected != actual:
    import sys, traceback
    (filename, lineno, container, code) = traceback.extract_stack()[-2]
    print("Test: %r failed on line %d in file %r.\nExpected %r but got %r\n" %
          (code, lineno, filename, expected, actual))

    sys.exit(1)


class TestEuclidean(unittest.TestCase):
  """"
  Basic tests for the extended euclidean algorithm
  """

  def test_gcd_ints(self):
    """Test that the GCD functions correctly on ints."""
    assert gcd(7, 9) == 1
    assert gcd(8, 18) == 2
    assert gcd(-12, 24) == -12
    # gcd is only unique up to multiplication by a unit, and so sometimes we'll get negatives.
    assert gcd(12, -24) == 12
    assert gcd(4864, 3458) == 38

  def test_extended_euclidean_ints(self):
    """Test that the extended Euclidean algorithm works correctly on ints."""
    assert extended_euclidean_algorithm(4864, 3458) == (32, -45, 38)
    assert 32 * 4864- 45 * 3458 == 38
    assert extended_euclidean_algorithm(3458, 4864) == (-45, 32, 38)

