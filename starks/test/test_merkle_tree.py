import unittest
from starks.modp import IntegersModP
from starks.merkle_tree import merkelize
from starks.merkle_tree import mk_branch 
from starks.merkle_tree import verify_branch

class TestMerkleTree(unittest.TestCase):
  """
  Basic tests for Merkle trees. 

  Mostly to build understanding of merkle-tree functioning. 
  """

  def test_merkletree(self):
    """Constructs and tests merkle tree"""
    t = merkelize([x.to_bytes(32, 'big') for x in range(128)])
    b = mk_branch(t, 59)
    # Is t[1] the merkle root? Shouldn't it be t[0]?
    assert verify_branch(t[1], 59, b, output_as_int=True) == 59
    print('Merkle tree works')

  def test_merkletree_zmodp(self):
    """Constructs merkle tree of field elements."""
    modulus = 7
    mod7 = IntegersModP(modulus)
    l = [mod7(i) for i in range(128)]
    m_tree = merkelize(l)
    assert len(m_tree) == 256

  def test_mk_branch(self):
    """Tests construction of a merkle tree branch."""
    t = merkelize([x.to_bytes(32, 'big') for x in range(128)])
    assert len(t) == 256
    b = mk_branch(t, 59)
    assert len(b) == 8

    t = merkelize([x.to_bytes(32, 'big') for x in range(256)])
    assert len(t) == 512
    b = mk_branch(t, 59)
    assert len(b) == 9

  def test_verify_branch(self):
    """Tests the verification of a merkle tree branch."""
    t = merkelize([x.to_bytes(32, 'big') for x in range(128)])
    assert len(t) == 256
    b = mk_branch(t, 59)
    assert len(b) == 8
    v = verify_branch(t[1], 59, b, output_as_int=True) == 59
    assert v
