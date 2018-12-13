try:
  from hashlib import blake2s
except:
  from pyblake2 import blake2s
blake = lambda x: blake2s(x).digest()


def permute4(values):
  """Permutes the list by selecting elements a quarter down 
  
  TODO(rbharath): The structure of this permutation seems
  rather arbitrary. Once I understand the codebase better,
  worth refactoring.
  """
  o = []
  ld4 = len(values) // 4
  for i in range(ld4):
    o.extend(
        [values[i], values[i + ld4], values[i + ld4 * 2], values[i + ld4 * 3]])
  return o


def get_index_in_permuted(x, L):
  """Mapping to undo the effects of permute4
  
  Useful to find the location of an index in original list
  in Merkle tree.
  """
  ld4 = L // 4
  return x // ld4 + 4 * (x % ld4)


def merkelize(L):
  """Creates a merkle-tree representation of the given list.
  
  The merkle-tree is stored as a list of length 2*len(L).
  The last len(L) elements are the original list elements
  (converted to 32 byte representations). The next L-1
  eleemnts store Merkle-tree elements. nodes[1] is the root
  of the merkle tree.
  """
  L = permute4(L)
  nodes = [b''] * len(L) + [
      x.to_bytes(32, 'big') if isinstance(x, int) else x for x in L
  ]
  for i in range(len(L) - 1, 0, -1):
    nodes[i] = blake(nodes[i * 2] + nodes[i * 2 + 1])
  return nodes


def mk_branch(tree, index):
  """A branch of the merkle tree is a list"""
  index = get_index_in_permuted(index, len(tree) // 2)
  index += len(tree) // 2
  o = [tree[index]]
  while index > 1:
    # I think this is to capture the sibling nodes.
    o.append(tree[index ^ 1])
    index //= 2
  return o


def verify_branch(root, index, proof, output_as_int=False):
  index = get_index_in_permuted(index, 2**len(proof) // 2)
  # I think this is a bug. Should be 2**len(proof) // 2
  # But I think it's OK since we're only doing parity checks
  #index += 2**len(proof)
  index += 2**len(proof) // 2
  v = proof[0]
  for p in proof[1:]:
    if index % 2:
      v = blake(p + v)
    else:
      v = blake(v + p)
    index //= 2
  assert v == root
  return int.from_bytes(proof[0], 'big') if output_as_int else proof[0]
