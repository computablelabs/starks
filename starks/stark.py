import time

modulus = 2**256 - 351 * 2**32 + 1

#class PrimeField():
#    """A simple class implementing a prime field."""
#
#    def __init__(self, modulus):
#        # Quick primality test
#        assert pow(2, modulus, modulus) == 2
#        self.modulus = modulus
#
#    def add(self, x, y):
#        return (x+y) % self.modulus
#
#    def sub(self, x, y):
#        return (x-y) % self.modulus
#
#    def mul(self, x, y):
#        return (x*y) % self.modulus
#
#    def inv(self, a):
#        """Modular inverse using the extended Euclidean algorithm"""
#        if a == 0:
#            return 0
#        lm, hm = 1, 0
#        low, high = a % self.modulus, self.modulus
#        while low > 1:
#            r = high//low
#            nm, new = hm-lm*r, high-low*r
#            lm, low, hm, high = nm, new, lm, low
#        return lm % self.modulus
#
#    def multi_inv(self, values):
#        partials = [1]
#        for i in range(len(values)):
#            partials.append(self.mul(partials[-1], values[i] or 1))
#        inv = self.inv(partials[-1])
#        outputs = [0] * len(values)
#        for i in range(len(values), 0, -1):
#            outputs[i-1] = self.mul(partials[i-1], inv) if values[i-1] else 0
#            inv = self.mul(inv, values[i-1] or 1)
#        return outputs
#
#    # Evaluate a polynomial at a point
#    def eval_poly_at(self, p, x):
#        y = 0
#        power_of_x = 1
#        for i, p_coeff in enumerate(p):
#            y += power_of_x * p_coeff
#            power_of_x = (power_of_x * x) % self.modulus
#        return y % self.modulus
#
#    # Build a polynomial that returns 0 at all specified xs
#    def zpoly(self, xs):
#        root = [1]
#        for x in xs:
#            root.insert(0, 0)
#            for j in range(len(root)-1):
#                root[j] -= root[j+1] * x
#        return [x % self.modulus for x in root]
#
#    def lagrange_interp(self, xs, ys):
#        # Generate master numerator polynomial, eg. (x - x1) * (x - x2) * ... * (x - xn)
#        root = self.zpoly(xs)
#
#        # Generate per-value numerator polynomials, eg. for x=x2,
#        # (x - x1) * (x - x3) * ... * (x - xn), by dividing the master
#        # polynomial back by each x coordinate
#        nums = [self.div_polys(root, [-x, 1]) for x in xs]
#
#        # Generate denominators by evaluating numerator polys at each x
#        denoms = [self.eval_poly_at(nums[i], xs[i]) for i in range(len(xs))]
#        invdenoms = self.multi_inv(denoms)
#
#        # Generate output polynomial, which is the sum of the per-value numerator
#        # polynomials rescaled to have the right y values
#        b = [0 for y in ys]
#        for i in range(len(xs)):
#            yslice = self.mul(ys[i], invdenoms[i])
#            for j in range(len(ys)):
#                if nums[i][j] and ys[i]:
#                    b[j] += nums[i][j] * yslice
#        return [x % self.modulus for x in b]

#def mimc(inp, steps, round_constants):
#    start_time = time.time()
#    for i in range(steps-1):
#        inp = (inp**3 + round_constants[i % len(round_constants)]) % modulus
#    print("MIMC computed in %.4f sec" % (time.time() - start_time))
#    return inp
