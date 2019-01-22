import unittest
from starks import blum_blum_shub
from starks.zkp_3_coloring import ZK3ColProver
from starks.zkp_3_coloring import ZK3ColVerifier

class TestZeroKnowledge3Coloring(unittest.TestCase):
  """
  Basic tests for zero knowledge 3 coloring statements.
  """

  def test_basic_example(self):
    """Basic tests"""
    # a graph is a list of edges, and for simplicity we'll say
    # every vertex shows up in some edge
    exampleGraph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    """
    A 3-coloring is a {int: int} where the output int is 0, 1, or 2.
    Note that we want to have as few bits as possible, since the bit
    commitment scheme blows up the size of the coloring by a factor
    of n.
    """
    exampleColoring = {
        1: 0,
        2: 1,
        3: 2,
        4: 1,
        5: 2,
        6: 0,
    }

  def test_commit_to_coloring(self):
    """Test commit to coloring"""
    security_parameter = 512
    G = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    coloring = {
        1: 0,
        2: 1,
        3: 2,
        4: 1,
        5: 2,
        6: 0,
    }
    one_way_permutation = blum_blum_shub.blum_blum_shub(security_parameter)
    hardcore_predicate = blum_blum_shub.parity
    prover = ZK3ColProver(G, coloring, one_way_permutation, hardcore_predicate)
    committed_coloring = prover.commit_to_coloring()
    assert len(committed_coloring) == 6

  def test_basic(self):
    """Conducts a basic test of the running the 3-coloring protocol."""
    security_parameter = 512
    one_way_permutation = blum_blum_shub.blum_blum_shub(security_parameter)
    hardcore_predicate = blum_blum_shub.parity
    G = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    coloring = {
        1: 0,
        2: 1,
        3: 2,
        4: 1,
        5: 2,
        6: 0,
    }

    prover = ZK3ColProver(G, coloring, one_way_permutation, hardcore_predicate)
    verifier = ZK3ColVerifier(G, one_way_permutation, hardcore_predicate)

    committed_coloring = prover.commit_to_coloring()
    chosen_edge = verifier.choose_edge(committed_coloring)

    revealed = prover.reveal_colors(*chosen_edge)
    revealed_colors = (
        verifier.verifier.decode(revealed[0], committed_coloring[chosen_edge[0]]),
        verifier.verifier.decode(revealed[1], committed_coloring[chosen_edge[1]]),
    )
    is_valid = verifier.accepts(revealed)

    print("{} != {} and commitment is valid? {}".format(
        revealed_colors[0], revealed_colors[1], is_valid))

