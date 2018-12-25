import random

from starks.zero_knowledge import num_vertices
from starks.zero_knowledge import random_permutation
from starks import blum_blum_shub
from starks import commitment

ONE_WAY_PERMUTATION = blum_blum_shub.blum_blum_shub(512)
HARDCORE_PREDICATE = blum_blum_shub.parity


class ZK3ColProver(object):
  """Proves in zero-knowledge a 3-coloring for a given graph.

  TODO(rbharath): This class is pretty stateful. Swap out for a less stateful
  implementation.
  """

  def __init__(self,
               graph,
               coloring,
               one_way_permutation=ONE_WAY_PERMUTATION,
               hardcore_predicate=HARDCORE_PREDICATE):
    self.graph = [tuple(sorted(e)) for e in graph]
    self.coloring = coloring
    self.vertices = list(range(1, num_vertices(graph) + 1))
    self.one_way_permutation = one_way_permutation
    self.hardcore_predicate = hardcore_predicate
    self.vertexToScheme = None

  def commit_to_coloring(self):
    """Commit to a coloring scheme"""

    self.vertexToScheme = {
        v: commitment.BBSIntCommitmentScheme(2, self.one_way_permutation,
                                             self.hardcore_predicate)
        for v in self.vertices
    }

    permutation = random_permutation(3)
    permutedColoring = {v: permutation[self.coloring[v]] for v in self.vertices}

    return {
        v: s.commit(permutedColoring[v])
        for (v, s) in self.vertexToScheme.items()
    }

  def reveal_colors(self, u, v):
    u, v = min(u, v), max(u, v)
    if not (u, v) in self.graph:
      raise Exception('Must query an edge!')

    return (
        self.vertexToScheme[u].reveal(),
        self.vertexToScheme[v].reveal(),
    )


class ZK3ColVerifier(object):
  """Verifies a zero-knowledge proof of 3-coloring."""

  def __init__(self,
               graph,
               one_way_permutation=ONE_WAY_PERMUTATION,
               hardcore_predicate=HARDCORE_PREDICATE):
    self.graph = [tuple(sorted(e)) for e in graph]
    self.one_way_permutation = one_way_permutation
    self.hardcore_predicate = hardcore_predicate
    self.committed_coloring = None
    self.verifier = commitment.BBSIntCommitmentVerifier(2, one_way_permutation,
                                                        hardcore_predicate)

  def choose_edge(self, committed_coloring):
    self.committed_coloring = committed_coloring
    self.chosen_edge = random.choice(self.graph)
    return self.chosen_edge

  def accepts(self, revealed):
    revealed_colors = []

    for (w, bitSecrets) in zip(self.chosen_edge, revealed):
      trueColor = self.verifier.decode(bitSecrets, self.committed_coloring[w])
      revealed_colors.append(trueColor)
      if not self.verifier.verify(bitSecrets, self.committed_coloring[w]):
        return False

    return revealed_colors[0] != revealed_colors[1]


def run_protocol(G, coloring, security_parameter=512):
  """Runs one round of the 3 coloring protocol for two graphs."""
  one_way_permutation = blum_blum_shub.blum_blum_shub(security_parameter)
  hardcore_predicate = blum_blum_shub.parity

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

  return is_valid
