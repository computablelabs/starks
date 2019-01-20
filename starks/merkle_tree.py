try:
  from hashlib import blake2s
except:
  from pyblake2 import blake2s
blake = lambda x: blake2s(x).digest()
from typing import List
from starks.numbertype import FieldElement
from starks.numbertype import Poly


def permute4(values: List) -> List:
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


def merkelize(L: List[FieldElement]) -> List[bytes]:
  """Creates a merkle-tree representation of the given list.
  
  The merkle-tree is stored as a list of length 2*len(L).
  The last len(L) elements are the original list elements
  (converted to 32 byte representations). The next L-1
  eleemnts store Merkle-tree elements. nodes[1] is the root
  of the merkle tree.
  """
  L = permute4(L)
  nodes = [b''] * len(L)
  for x in L:
    if isinstance(x, int):
      nodes.append(x.to_bytes(32, 'big'))
    elif isinstance(x, bytes):
      nodes.append(x)
    else:
      nodes.append(x.to_bytes())
  for i in range(len(L) - 1, 0, -1):
    nodes[i] = blake(nodes[i * 2] + nodes[i * 2 + 1])
  return nodes


def mk_branch(tree: List[FieldElement], index: int) -> List[FieldElement]:
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
  """Verifies the proof and returns the leaf on the branch"""
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

def evaluate_polynomials(polynomials: List[Poly]):
  """Convert polynomials into evaluated form."""
  fft = MultiDimNonBinaryFFT(field, root_of_unity, width)
  values = fft.multi_fft(polynomials)
  return values

def merkelize_polynomial_evaluations(dims, polynomial_evals: List[List[FieldElement]]):
  """Given a list of polynomial evaluations, merkelizes them together.

  Each leaf of the Merkle tree contains the concatenation of the values of
  each polynomial in polynomials at a given index. If the polynomials are
  multidimensional, the per-dimension leaves are concatenated to form one
  joint leaf.

  Parameters
  ----------
  dims: Int
    Dimensionality
  polynomials: List of polynomials by dimension
    Each element much be a list of evaluations of a given poly. All of
    these should have the same length.
  """
  # Note len(mtree_leaf) == 32 * dims * len(polys)
  # Note len(mtree) == 2 * precision now
  # This code packs each Merkle leaf as
  # polyval_1_dim_1 polyval_1_dim_2 poly_val_2_dim_1 poly_val_2_dim_2 ...
  # In the common case this is
  # [p_of_x_dim_1 p_of_x_dim_2 .. d_of_x_dim_1 d_of_x_dim_2... b_of_x_dim_1 b_of_x_dim_2]
  ###################################################
  #print("type(polynomial_evals[0][0])")
  #print(type(polynomial_evals[0][0]))
  #print("polynomial_evals[0][0].to_bytes()")
  #print(polynomial_evals[0][0].to_bytes())
  #for evals in zip(*polynomial_evals):
  #  print("evals")
  #  print(evals)
  #  for val in evals:
  #    print("type(val)")
  #    print(type(val))
  #    print("val.to_bytes()")
  #    print(val.to_bytes())
  #  break
  ###################################################
  mtree = merkelize([
      # TODO(rbharath): Assuming now in field. May fix later
      #b''.join([val[dim].to_bytes(32, 'big') for val in evals for dim in range(dims)])
      #b''.join([val[dim].to_bytes() for val in evals for dim in range(dims)])
      b''.join([val.to_bytes() for val in evals])
      for evals in zip(*polynomial_evals)])
  return mtree

def unpack_merkle_leaf(leaf: bytes, dims: int, num_polys: int) -> List[bytes]:
  """Unpacks a merkle leaf created by merkelize_polynomials.

  Note that the packing in each Merkle leaf is

  polyval_1_dim_1 polyval_1_dim_2 poly_val_2_dim_1 poly_val_2_dim_2 ...

  Each value is encoded in 32 bytes, so the total length of
  the leaf is 32 * num_polys * dims bytes.

  Parameters
  ----------
  leaf: bytestring
    A bytestring holding the Merkle leaf
  dims: Int
    The dimensionality of the state space
  num_polys: int
    The number of polynomials
  """
  vals = []
  for poly_ind in range(num_polys):
    for dim in range(dims):
      start_index = 32 * (poly_ind * dims + dim)
      end_index = 32 * (poly_ind * dims + dim + 1)
      byte_val = leaf[start_index:end_index]
      vals.append(byte_val)
  return vals
