import pytest
from starks.merkle_tree import merkelize
from starks.merkle_tree import mk_branch 
from starks.merkle_tree import verify_branch
#from starks.modp import IntegersModP
#from starks.rationals_modp import RationalsModP
#from starks.merkle_tree import unpack_merkle_leaf

def test_merkletree():
    """Constructs and tests merkle tree"""
    t = merkelize([x.to_bytes(32, 'big') for x in range(128)])
    b = mk_branch(t, 59)
    # Is t[1] the merkle root? Shouldn't it be t[0]?
    assert verify_branch(t[1], 59, b, output_as_int=True) == 59
    print('Merkle tree works')
#
#def test_merkletree_zmodp():
#  """Constructs merkle tree of field elements."""
#  modulus = 7
#  mod7 = IntegersModP(modulus)
#  l = [mod7(i) for i in range(128)]
#  m_tree = merkelize(l)
#  assert len(m_tree) == 256
#
#def test_merkletree_qmodp():
#  """Constructs merkle tree of field elements."""
#  modulus = 7
#  ratmod = RationalsModP(modulus)
#  l = [ratmod(i, j) for i in range(12) for j in range(12)]
#  m_tree = merkelize(l)
#  assert len(m_tree) == 288
#
#
#def test_mk_branch():
#  """Tests construction of a merkle tree branch."""
#  t = merkelize([x.to_bytes(32, 'big') for x in range(128)])
#  assert len(t) == 256
#  b = mk_branch(t, 59)
#  assert len(b) == 8
#
#  t = merkelize([x.to_bytes(32, 'big') for x in range(256)])
#  assert len(t) == 512
#  b = mk_branch(t, 59)
#  assert len(b) == 9
#
#def test_verify_branch():
#  """Tests the verification of a merkle tree branch."""
#  t = merkelize([x.to_bytes(32, 'big') for x in range(128)])
#  assert len(t) == 256
#  b = mk_branch(t, 59)
#  assert len(b) == 8
#  v = verify_branch(t[1], 59, b, output_as_int=True) == 59
#  assert v
#
#def test_unpack_merkle_leaf():
#  """
#  Tests that merkle leaf unpacking works correctly.
#  """
#  leaf_parts = []
#  dims = 2
#  num_polys = 3
#  count = 0
#  for poly_ind in range(num_polys):
#    for dim in range(dims):
#      count_bytes = count.to_bytes(32, 'big')
#      leaf_parts.append(count_bytes)
#      count += 1
#  leaf = b''.join(leaf_parts)
#  vals = unpack_merkle_leaf(leaf, dims, num_polys)
#  assert len(vals) == len(leaf_parts)
#  for ind, val in enumerate(vals):
#    assert val == leaf_parts[ind]
#
