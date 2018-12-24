import random

from starks.zero_knowledge import num_vertices
from starks.zero_knowledge import random_permutation
import blum_blum_shub
import commitment

ONE_WAY_PERMUTATION = blum_blum_shub.blum_blum_shub(512)
HARDCORE_PREDICATE = blum_blum_shub.parity


class Prover(object):

  def __init__(self,
               graph,
               coloring,
               oneWayPermutation=ONE_WAY_PERMUTATION,
               hardcorePredicate=HARDCORE_PREDICATE):
    self.graph = [tuple(sorted(e)) for e in graph]
    self.coloring = coloring
    self.vertices = list(range(1, num_vertices(graph) + 1))
    self.oneWayPermutation = oneWayPermutation
    self.hardcorePredicate = hardcorePredicate
    self.vertexToScheme = None

  def commitToColoring(self):
    self.vertexToScheme = {
        v: commitment.BBSIntCommitmentScheme(2, self.oneWayPermutation,
                                             self.hardcorePredicate)
        for v in self.vertices
    }

    permutation = random_permutation(3)
    permutedColoring = {v: permutation[self.coloring[v]] for v in self.vertices}

    return {
        v: s.commit(permutedColoring[v])
        for (v, s) in self.vertexToScheme.items()
    }

  def revealColors(self, u, v):
    u, v = min(u, v), max(u, v)
    if not (u, v) in self.graph:
      raise Exception('Must query an edge!')

    return (
        self.vertexToScheme[u].reveal(),
        self.vertexToScheme[v].reveal(),
    )


class Verifier(object):

  def __init__(self,
               graph,
               oneWayPermutation=ONE_WAY_PERMUTATION,
               hardcorePredicate=HARDCORE_PREDICATE):
    self.graph = [tuple(sorted(e)) for e in graph]
    self.oneWayPermutation = oneWayPermutation
    self.hardcorePredicate = hardcorePredicate
    self.committedColoring = None
    self.verifier = commitment.BBSIntCommitmentVerifier(2, oneWayPermutation,
                                                        hardcorePredicate)

  def chooseEdge(self, committedColoring):
    self.committedColoring = committedColoring
    self.chosenEdge = random.choice(self.graph)
    return self.chosenEdge

  def accepts(self, revealed):
    revealedColors = []

    for (w, bitSecrets) in zip(self.chosenEdge, revealed):
      trueColor = self.verifier.decode(bitSecrets, self.committedColoring[w])
      revealedColors.append(trueColor)
      if not self.verifier.verify(bitSecrets, self.committedColoring[w]):
        return False

    return revealedColors[0] != revealedColors[1]


def runProtocol(G, coloring, securityParameter=512):
  oneWayPermutation = blum_blum_shub.blum_blum_shub(securityParameter)
  hardcorePredicate = blum_blum_shub.parity

  prover = Prover(G, coloring, oneWayPermutation, hardcorePredicate)
  verifier = Verifier(G, oneWayPermutation, hardcorePredicate)

  committedColoring = prover.commitToColoring()
  chosenEdge = verifier.chooseEdge(committedColoring)

  revealed = prover.revealColors(*chosenEdge)
  revealedColors = (
      verifier.verifier.decode(revealed[0], committedColoring[chosenEdge[0]]),
      verifier.verifier.decode(revealed[1], committedColoring[chosenEdge[1]]),
  )
  isValid = verifier.accepts(revealed)

  print("{} != {} and commitment is valid? {}".format(
      revealedColors[0], revealedColors[1], isValid))

  return isValid


if __name__ == "__main__":
  for _ in range(30):
    runProtocol(exampleGraph, exampleColoring, securityParameter=10)
