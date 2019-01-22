"""Basic tests for commitment scheme."""
import unittest
from starks import blum_blum_shub
from starks.commitment import BBSBitCommitmentScheme
from starks.commitment import BBSBitCommitmentVerifier

class TestCommitment(unittest.TestCase):
  """"
  Basic tests for Commitment schemes. 
  """

  def test_pred_construction(self):
    """Test predicate construction works."""
    security_parameter = 10
    one_way_perm = blum_blum_shub.blum_blum_shub(security_parameter)
    hardcorePred = blum_blum_shub.parity

  def test_bit_commitment(self):
    """Test that bit commitment works."""
    security_parameter = 10
    one_way_perm = blum_blum_shub.blum_blum_shub(security_parameter)
    hardcorePred = blum_blum_shub.parity

    print('Bit commitment')
    scheme = BBSBitCommitmentScheme(one_way_perm, hardcorePred, security_parameter)
    verifier = BBSBitCommitmentVerifier(one_way_perm, hardcorePred)
