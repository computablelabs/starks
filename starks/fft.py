def _simple_ft(vals, modulus, roots_of_unity):
  """Efficient base case implementation.
  
  The FFT recurses down halves of the list. This method is
  called to handle the base case of the fft.
  """
  L = len(roots_of_unity)
  o = []
  for i in range(L):
    last = 0
    for j in range(L):
      last += vals[j] * roots_of_unity[(i * j) % L]
    o.append(last % modulus)
  return o


def _fft(vals, modulus, roots_of_unity):
  if len(vals) <= 4:
    #return vals
    return _simple_ft(vals, modulus, roots_of_unity)
  L = _fft(vals[::2], modulus, roots_of_unity[::2])
  R = _fft(vals[1::2], modulus, roots_of_unity[::2])
  o = [0 for i in vals]
  for i, (x, y) in enumerate(zip(L, R)):
    y_times_root = y * roots_of_unity[i]
    o[i] = (x + y_times_root) % modulus
    o[i + len(L)] = (x - y_times_root) % modulus
  return o

def fft(vals, modulus, root_of_unity, inv=False, dims=1):
  """Computes FFT for potentially multidimensional sequences"""
  if dims == 1:
    return fft_1d(vals, modulus, root_of_unity, inv)
  else:
    fft_vals = []
    for dim in range(dims):
      vals_dim = [val[dim] for val in vals]
      fft_dim = fft_1d(vals_dim, modulus, root_of_unity, inv=inv)
      fft_vals.append(fft_dim)
    fft_joint = list(zip(*fft_vals))
    return fft_joint


def fft_1d(vals, modulus, root_of_unity, inv=False):
  """Computes FFT for one dimensional inputs"""
  # Build up roots of unity
  rootz = [1, root_of_unity]
  while rootz[-1] != 1:
    rootz.append((rootz[-1] * root_of_unity) % modulus)
  # Fill in vals with zeroes if needed
  if len(rootz) > len(vals) + 1:
    vals = vals + [0] * (len(rootz) - len(vals) - 1)
  if inv:
    # Inverse FFT
    invlen = pow(len(vals), modulus - 2, modulus)
    return [(x * invlen) % modulus for x in _fft(vals, modulus, rootz[:0:-1])]
  else:
    # Regular FFT
    return _fft(vals, modulus, rootz[:-1])


def mul_polys(a, b, modulus, root_of_unity):
  """Multiply polynomials by converting to fourier space"""
  rootz = [1, root_of_unity]
  while rootz[-1] != 1:
    rootz.append((rootz[-1] * root_of_unity) % modulus)
  if len(rootz) > len(a) + 1:
    a = a + [0] * (len(rootz) - len(a) - 1)
  if len(rootz) > len(b) + 1:
    b = b + [0] * (len(rootz) - len(b) - 1)
  x1 = _fft(a, modulus, rootz[:-1])
  x2 = _fft(b, modulus, rootz[:-1])
  return _fft([(v1 * v2) % modulus for v1, v2 in zip(x1, x2)], modulus,
              rootz[:0:-1])
