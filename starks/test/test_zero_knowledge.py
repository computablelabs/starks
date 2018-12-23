import unittest
import random
import types
from starks.zero_knowledge import num_vertices
from starks.zero_knowledge import random_permutation
from starks.zero_knowledge import make_permutation_function
from starks.zero_knowledge import make_inverse_permutation_function
from starks.zero_knowledge import apply_isomorphism
from starks.zero_knowledge import run_protocol
from starks.zero_knowledge import convince_beyond_doubt
from starks.zero_knowledge import ZKProver
from starks.zero_knowledge import ZKVerifier

class TestZeroKnowledge(unittest.TestCase):
  """
  Basic tests for zero knowledge handling.
  """

  def test_graph_ops(self):
    """
    Basic tests for graph handling.
    """
    # a graph is a list of edges, and for simplicity we'll say
    # every vertex shows up in some edge
    exampleGraph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    assert num_vertices(exampleGraph) == 6

  def test_permutation_ops(self):
    """
    Basic tests for how permutations work.
    """
    example_graph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    perm_list = random_permutation(6)
    perm_f = make_permutation_function(perm_list)
    permuted_graph = apply_isomorphism(example_graph, perm_f)
    assert len(permuted_graph) == len(example_graph)
    assert num_vertices(permuted_graph) == num_vertices(example_graph)

  def test_zk_graph_prover_construction(self):
    """
    Basic tests on how the ZK graph prover works.
    """
    example_graph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    perm_list = random_permutation(6)
    perm_f = make_permutation_function(perm_list)
    permuted_graph = apply_isomorphism(example_graph, perm_f)
    prover = ZKProver(example_graph, permuted_graph, perm_list)

  def test_zk_graph_prover_send_isomorphic(self):
    """
    Testing generation of random isomorphic graph. 
    """
    example_graph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    perm_list = random_permutation(6)
    perm_f = make_permutation_function(perm_list)
    permuted_graph = apply_isomorphism(example_graph, perm_f)
    prover = ZKProver(example_graph, permuted_graph, perm_list)
    H = prover.send_isomorphic_copy()
    assert num_vertices(H) == num_vertices(example_graph)
    assert len(H) == len(example_graph)

  def test_zk_graph_prove_isomorphic_to(self):
    """
    Test proofs of graph isomorphism.

    TODO(rbharath): Create a better test for this function.
    """
    example_graph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    perm_list = random_permutation(6)
    perm_f = make_permutation_function(perm_list)
    permuted_graph = apply_isomorphism(example_graph, perm_f)
    prover = ZKProver(example_graph, permuted_graph, perm_list)
    H = prover.send_isomorphic_copy()
    iso_1 = prover.prove_isomorphic_to(1)
    assert isinstance(iso_1, types.FunctionType)
    iso_2 = prover.prove_isomorphic_to(2)
    assert isinstance(iso_2, types.FunctionType)

  def test_zk_verifier_construction(self):
    """Tests construction of verifier."""
    example_graph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    perm_list = random_permutation(6)
    perm_f = make_permutation_function(perm_list)
    permuted_graph = apply_isomorphism(example_graph, perm_f)
    verifier = ZKVerifier(example_graph, permuted_graph)

  def test_zk_verifier_choose_graph(self):
    """Tests verifier graph choice."""
    example_graph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    perm_list = random_permutation(6)
    perm_f = make_permutation_function(perm_list)
    permuted_graph = apply_isomorphism(example_graph, perm_f)
    prover = ZKProver(example_graph, permuted_graph, perm_list)
    H = prover.send_isomorphic_copy()
    verifier = ZKVerifier(example_graph, permuted_graph)
    choice = verifier.choose_graph(H)
    assert choice in [1, 2]

  def test_zk_verifier_acceptance(self):
    """Tests verifier graph accepts correct."""
    example_graph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    perm_list = random_permutation(6)
    perm_f = make_permutation_function(perm_list)
    permuted_graph = apply_isomorphism(example_graph, perm_f)
    prover = ZKProver(example_graph, permuted_graph, perm_list)
    H = prover.send_isomorphic_copy()
    verifier = ZKVerifier(example_graph, permuted_graph)
    choice = verifier.choose_graph(H)
    witness = prover.prove_isomorphic_to(choice)
    assert verifier.accepts(witness)

  # TODO(rbharath): This breaks!! Why??
  def test_zk_verifier_id_acceptance(self):
    G1 = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    G2 = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6), (1,6)]
    perm_list = [1, 2, 3, 4, 5, 6]
    prover = ZKProver(G1, G2, perm_list)
    verifier = ZKVerifier(G1, G2)

    H = prover.send_isomorphic_copy()
    choice = verifier.choose_graph(H)
    witness = prover.prove_isomorphic_to(choice)
    assert verifier.accepts(witness)

  def test_zk_verifier_nonacceptance(self):
    """Tests verifier rejects wrong isomorphism."""
    # Two non-isomorphic graphs
    G1 = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    G2 = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6), (1,6)]
    perm_list = [1, 2, 3, 4, 5, 6]
    trials = 10
    results = []
    success = False
    for i in range(trials):
      prover = ZKProver(G1, G2, perm_list)
      H = prover.send_isomorphic_copy()
      verifier = ZKVerifier(G1, G2)
      choice = verifier.choose_graph(H)
      witness = prover.prove_isomorphic_to(choice)
      result = verifier.accepts(witness)
      if not result:
        success = True
      results.append(result)
    assert success

  def test_run_protocol(self):
    """Test a run of the protocol."""
    G1 = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    perm_list = random_permutation(6)
    perm_f = make_permutation_function(perm_list)
    G2 = apply_isomorphism(G1, perm_f)
    result = run_protocol(G1, G2, perm_list)
    print("result")
    print(result)

  def test_basic(self):
    """Basic test in original source."""
    example_graph = [(1, 2), (1, 4), (1, 3), (2, 5), (2, 5), (3, 6), (5, 6)]
    G1 = example_graph
    n = num_vertices(G1)
    p = random_permutation(n)

    f = make_permutation_function(p)
    finv = make_inverse_permutation_function(p)
    G2 = apply_isomorphism(G1, f)

    assert apply_isomorphism(G1, f) == G2
    assert apply_isomorphism(G2, finv) == G1
