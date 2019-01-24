"""Tests for the Blum-blum-shub pseudorandom generator."""

import unittest
from starks.blum_blum_shub import blum_blum_shub
from starks.blum_blum_shub import parity

class TestBBS(unittest.TestCase):
  """"
  Basic tests for Blum-blum-shub 
  """

  # TODO(rbharath): Expand these tests.
  def test_basic(self):
    """Basic tests."""
    owp = blum_blum_shub()
    print(owp(70203203))
    print(owp(12389))

  def test_parity(self):
    """Tests the parity operation."""
    # parity(3) == 0 since 3 = 0b11 has 2 one's.
    assert parity(3) == 0
    # parity(4) == 1 since 4 = 0b100 has 1 one.
    assert parity(4) == 1


