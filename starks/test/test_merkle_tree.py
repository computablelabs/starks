import unittest
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

