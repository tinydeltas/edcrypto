""""
Base class for all modular arithmetic operations.
"""

import random
import ecdsa.numbertheory
from group import Secret

class ModInt(Secret, object):
    def __init__(self, p=None, v=None):
        """
        p - prime
        v - actual value of the secret
        """
        self.p = p
        self.v = v

    def string(self):
        return self.v

    def zero(self):
        self.v = 0
        return self

    def one(self):
        self.v = 1
        return self

    def add(self, a, b):
        self.v = (a.v + b.v) % self.p.v
        return self

    def sub(self, a, b):
        self.v = (a.v - b.v) % self.p.v
        return self

    def neg(self, a):
        self.v = -a.v % self.p.v
        return self

    def mul(self, a, b):
        self.v = (a.v * b.v) % self.p.v
        return self

    def equal(self, b):
        """
        Takes two modular ints and checks them for equality
        """
        return (self.v % self.p.v) == (b.v % self.p.v)

    def set(self, a):
        """
        Sets the value of the int to a's value.
        """
        self.v = a.v
        self.p = a.p
        return self

    def inv(self, a):
        self.v = ecdsa.numbertheory.inverse_mod(a.v, self.p.v)
        return self

    def sqrt(self, a):
        self.v = ecdsa.numbertheory.square_root_mod_prime(a.v, self.p.v)
        return self

    def exp(self, a, exponent):
        self.v = ecdsa.numbertheory.modular_exp(a.v, exponent.v, self.p.v)
        return self

    def jacobi(self, a):
        return  ecdsa.numbertheory.jacobi(a.v, self.p.v)

    def div(self, a, b):
        self.v = (a.v * ecdsa.numbertheory.inverse_mod(b.v, self.p.v)) % self.p.v
        return self

    def random_secret(self):
        """
        Generates random secret in the range [1, self.p - 1]
        """
        self.v = random.randrange(1, self.p.v)
        return self

