import random

def num_vertices(G):
  """Counts the number of vertices in the graph."""
  return max(v for e in G for v in e)


def random_permutation(n):
  """Constructs a random permutation for n elements"""
  L = list(range(n))
  random.shuffle(L)
  return L


def make_permutation_function(L):
  """Transforms a list representation of permutation into function"""
  return lambda i: L[i - 1] + 1


def make_inverse_permutation_function(L):
  """Constructs inverse permutation function for given list"""
  return lambda i: 1 + L.index(i - 1)


def apply_isomorphism(G, f):
  """Applies the specified isomorphism to the graph."""
  return [(f(i), f(j)) for (i, j) in G]


class ZKProver(object):
  """Proves statements in zero-knowledge.

  At present, this class only works for proving graph isomorphism.s

  TODO(rbharath): Figure out how to extend this to computational
  transcripts. Applying an isomorphism to the trace doesn't quite work
  since we would need to encode the isomorphism in the step transition
  somehow.
  """

  def __init__(self, G1, G2, isomorphism):
    """
    Isomomorphism is a list of integers representing
    an isomorphism from G1 to G2.

    Parameters
    ----------
    G1: Graph
      First graph
    G2: Graph
      Second Graph
    isomorphism: function
      A function implementing an isomorphism
    """
    self.G1 = G1
    self.G2 = G2
    self.n = num_vertices(G1)
    assert self.n == num_vertices(G2)

    self.isomorphism = isomorphism
    self.state = None

  def send_isomorphic_copy(self):
    """Returns a permuted version of G1."""
    isomorphism = random_permutation(self.n)
    pi = make_permutation_function(isomorphism)

    H = apply_isomorphism(self.G1, pi)

    self.state = isomorphism
    return H

  def prove_isomorphic_to(self, graph_choice):
    """Prove isomorphic to either G1 or G2

    Must be called after a call to send_isomorphic_copy

    Parameters
    ----------
    graph_choice: int
      Either 1 or 2 for G1 and G2 respectively.
    """
    random_isomorphism = self.state
    pi_inverse = make_inverse_permutation_function(random_isomorphism)

    if graph_choice == 1:
      return pi_inverse
    else:
      f = make_permutation_function(self.isomorphism)
      return lambda i: f(pi_inverse(i))


class ZKVerifier(object):
  """Verifies that the prover is functioning correctly.

  TODO(rbharath): Create a better docstring.
  """

  def __init__(self, G1, G2):
    """
    Parameters
    ----------
    G1: graph
      First graph
    G2: graph
      Second graph
    """
    self.G1 = G1
    self.G2 = G2
    self.n = num_vertices(G1)
    assert self.n == num_vertices(G2)

  def choose_graph(self, H):
    """Makes a random choice between G1 and G2."""
    choice = random.choice([1, 2])
    self.state = H, choice
    return choice

  def accepts(self, isomorphism):
    """
    Return True if and only if the given isomorphism
    is a valid isomorphism between the randomly
    chosen graph in the first step, and the H presented
    by the Prover.
    """
    H, choice = self.state
    graph_to_check = [self.G1, self.G2][choice - 1]
    f = isomorphism

    # The correctness of this check depends on python list comparison being
    # done with == since each graph is a list of tuples.
    is_valid_isomorphism = (graph_to_check == apply_isomorphism(H, f))
    return is_valid_isomorphism


def run_protocol(G1, G2, isomorphism):
  """Runs one round of the protocol for two graphs."""
  p = ZKProver(G1, G2, isomorphism)
  v = ZKVerifier(G1, G2)

  H = p.send_isomorphic_copy()
  choice = v.choose_graph(H)
  witness_isomorphism = p.prove_isomorphic_to(choice)

  return v.accepts(witness_isomorphism)


def convince_beyond_doubt(G1, G2, isomorphism, error_tolerance=1e-20):
  """Runs the protocol sufficient times to prove without doubt."""
  probability_fooled = 1

  while probability_fooled > error_tolerance:
    result = run_protocol(G1, G2, isomorphism)
    assert result
    probability_fooled *= 0.5
    print(probability_fooled)


def messagesFromProtocol(G1, G2, isomorphism):
  p = Prover(G1, G2, isomorphism)
  v = Verifier(G1, G2)

  H = p.send_isomorphic_copy()
  choice = v.choose_graph(H)
  witnessIsomorphism = p.prove_isomorphic_to(choice)

  return [H, choice, witnessIsomorphism]


def simulateProtocol(G1, G2):
  # Construct data drawn from the same distribution as what is
  # returned by messagesFromProtocol
  choice = random.choice([1, 2])
  G = [G1, G2][choice - 1]
  n = num_vertices(G)

  isomorphism = random_permutation(n)
  pi = make_permutation_function(isomorphism)
  H = apply_isomorphism(G, pi)

  return H, choice, pi


if __name__ == "__main__":
  G1 = exampleGraph
  n = num_vertices(G1)
  p = random_permutation(n)

  f = make_permutation_function(p)
  finv = make_inverse_permutation_function(p)
  G2 = apply_isomorphism(G1, f)

  assert apply_isomorphism(G1, f) == G2
  assert apply_isomorphism(G2, finv) == G1

  convinceBeyondDoubt(G1, G2, p)
