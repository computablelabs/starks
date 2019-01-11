from typing import List
from starks.merkle_tree import merkelize
from starks.merkle_tree import mk_branch
from starks.merkle_tree import verify_branch
from starks.utils import get_power_cycle
from starks.utils import get_pseudorandom_indices
from starks.poly_utils import lagrange_interp
from starks.poly_utils import multi_interp_4
from starks.numbertype import FieldElement
from starks.numbertype import Field

# The number of spot checks performed at each recursion of the
# FRI proof.

def prove_low_degree(values: List[FieldElement],
                     root_of_unity: FieldElement,
                     maxdeg_plus_1: int,
                     #modulus: int,
                     field: Field,
                     exclude_multiples_of:int = 0,
                     fri_spot_check_security_factor:int = 40):
  """
  Generate an FRI proof that the polynomial that has the
  specified values at successive powers of the specified root
  of unity has a degree lower than maxdeg_plus_1
  
  We use maxdeg+1 instead of maxdeg because it's more
  mathematically convenient in this case.

  Note that if values is a n-degree polynomial, root_of_unity
  should be a n-th root of unity.
  """
  #f = PrimeField(modulus)
  print('Proving %d values are degree <= %d' % (len(values), maxdeg_plus_1))

  # If the degree we are checking for is less than or equal to
  # 32, use the polynomial directly as a proof
  # TODO(rbharath): Why does this make sense?
  if maxdeg_plus_1 <= 16:
    print('Produced FRI proof')
    #return [[x.to_bytes(32, 'big') for x in values]]
    return [[x.to_bytes() for x in values]]

  # Calculate the set of x coordinates
  #xs = get_power_cycle(root_of_unity, modulus)
  xs = get_power_cycle(root_of_unity, field)
  
  assert len(values) == len(xs)

  # Put the values into a Merkle tree. This is the root that
  # the proof will be checked against. Note this is a list of
  # length 2*len(values) storing the complete merkle tree.
  m = merkelize(values)

  # Select a pseudo-random x coordinate
  # This is the merkle-root of the polynomial.
  #special_x = int.from_bytes(m[1], 'big') % modulus
  special_x = field(m[1])

  # Calculate the "column" at that x coordinate (see
  # https://vitalik.ca/general/2017/11/22/starks_part_2.html)
  # We calculate the column by Lagrange-interpolating each
  # row, and not directly from the polynomial, as this is more
  # efficient
  quarter_len = len(xs) // 4
  x_polys = multi_interp_4(field,
      [[xs[i + quarter_len * j] for j in range(4)] for i in range(quarter_len)],
      [[values[i + quarter_len * j]
        for j in range(4)]
       for i in range(quarter_len)])
  #column = [f.eval_quartic(p, special_x) for p in x_polys]
  column = [p(special_x) for p in x_polys]
  m2 = merkelize(column)

  # Pseudo-randomly select y indices to sample
  ys = get_pseudorandom_indices(
      m2[1], len(column), fri_spot_check_security_factor, exclude_multiples_of=exclude_multiples_of)

  # Compute the Merkle branches for the values in the
  # polynomial and the column
  branches = []
  for y in ys:
    branches.append([mk_branch(m2, y)] +
                    [mk_branch(m, y + (len(xs) // 4) * j) for j in range(4)])

  # This component of the proof
  o = [m2[1], branches]

  # Recurse...
  return [o] + prove_low_degree(
      column,
      #f.exp(root_of_unity, 4),
      root_of_unity**4,
      maxdeg_plus_1 // 4,
      #modulus,
      field,
      exclude_multiples_of=exclude_multiples_of)


def verify_low_degree_proof(merkle_root,
                            root_of_unity,
                            proof,
                            maxdeg_plus_1,
                            #modulus,
                            field,
                            exclude_multiples_of=0,
                            fri_spot_check_security_factor=40):
  """Verify an FRI proof"""
  # Calculate which root of unity we're working with
  testval = root_of_unity
  # roudeg is the power of the root of unity 
  roudeg = 1
  while testval != 1:
    roudeg *= 2
    #testval = (testval * testval) % modulus
    testval = (testval * testval)

  # Powers of the given root of unity 1, p, p**2, p**3 such that p**4 = 1
  quartic_roots_of_unity = [
      1,
      root_of_unity**(roudeg // 4),
      root_of_unity**(roudeg // 2),
      root_of_unity**(roudeg * 3 // 4)
  ]

  # Verify the recursive components of the proof
  for prf in proof[:-1]:
    root2, branches = prf
    print('Verifying degree <= %d' % maxdeg_plus_1)

    # Calculate the pseudo-random x coordinate
    #special_x = int.from_bytes(merkle_root, 'big') % modulus
    special_x = field(merkle_root)

    # Calculate the pseudo-randomly sampled y indices
    ys = get_pseudorandom_indices(
        root2, roudeg // 4, fri_spot_check_security_factor, exclude_multiples_of=exclude_multiples_of)

    # For each y coordinate, get the x coordinates on the row,
    # the values on the row, and the value at that y from the
    # column
    xcoords = []
    rows = []
    columnvals = []
    for i, y in enumerate(ys):
      # The x coordinates from the polynomial
      x1 = root_of_unity**y
      xcoords.append(
          [(quartic_roots_of_unity[j] * x1) for j in range(4)])

      # The values from the original polynomial
      row = [
          verify_branch(
              merkle_root, y + (roudeg // 4) * j, prf, output_as_int=True)
          for j, prf in zip(range(4), branches[i][1:])
      ]
      rows.append(row)

      columnvals.append(
          verify_branch(root2, y, branches[i][0], output_as_int=True))

    # Verify for each selected y coordinate that the four
    # points from the polynomial and the one point from the
    # column that are on that y coordinate are on the same deg
    # < 4 polynomial
    #polys = multi_interp_4(modulus, xcoords, rows)
    polys = multi_interp_4(field, xcoords, rows)

    for p, c in zip(polys, columnvals):
      assert p(special_x) == c

    # Update constants to check the next proof
    merkle_root = root2
    root_of_unity = root_of_unity**4
    maxdeg_plus_1 //= 4
    roudeg //= 4

  # Verify the direct components of the proof
  data = [int.from_bytes(x, 'big') for x in proof[-1]]
  print('Verifying degree <= %d' % maxdeg_plus_1)
  assert maxdeg_plus_1 <= 16

  # Check the Merkle root matches up
  mtree = merkelize(data)
  assert mtree[1] == merkle_root

  # Check the degree of the data
  #powers = get_power_cycle(root_of_unity, modulus)
  powers = get_power_cycle(root_of_unity, field)
  if exclude_multiples_of:
    pts = [x for x in range(len(data)) if x % exclude_multiples_of]
  else:
    pts = range(len(data))

  #poly = f.lagrange_interp([powers[x] for x in pts[:maxdeg_plus_1]],
  #poly = lagrange_interp(modulus,
  poly = lagrange_interp(field,
          [powers[x] for x in pts[:maxdeg_plus_1]],
          [data[x] for x in pts[:maxdeg_plus_1]])
  for x in pts[maxdeg_plus_1:]:
    assert poly(powers[x]) == data[x]

  print('FRI proof verified')
  return True
