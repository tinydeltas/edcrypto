"""
Weierstrauss curve base class.
"""

import random
import math
import unittest
from numth import jacobi, sqrt, inverse
from utils import string_to_long, b2l, l2b

from elgamal import ElGamal, PublicKey, PrivateKey

# parameter initialization for NIST-endorsed p256 (prime field)
def curvep256():
    p = string_to_long("ffffffff 00000001 00000000 00000000 00000000 ffffffff"
        "ffffffff ffffffff")
    a = p - 3
    b = string_to_long("5ac635d8 aa3a93e7 b3ebbd55 769886bc 651d06b0 cc53b0f6"
        "3bce3c3e 27d2604b")
    q = string_to_long("ffffffff 00000000 ffffffff ffffffff bce6faad a7179e84"
        "f3b9cac2 fc632551")
    k = 1000
    g_x = string_to_long("6b17d1f2 e12c4247 f8bce6e5 63a440f2 77037d81"
        "2deb33a0 f4a13945 d898c296")
    g_y = string_to_long("4fe342e2 fe1a7f9b 8ee7eb4a 7c0f9e16 2bce3357"
        "6b315ece cbb64068 37bf51f5")
    return EllipticCurve(a, b, p, q, k, g_x, g_y)

# P-256 is represented in Weierstrauss formm, which has separate laws for
# addition and doubling.

class EllipticCurve(ElGamal):
    def __init__(self, a, b, p, q, k=None, g_x=None, g_y=None):
        self.a = a
        self.b = b
        self.p = p
        self.q = q
        self.k = k
        self.g = (g_x, g_y)

    def zero(self):
        return 0

    def generator(self):
        return self.g

    def order(self):
        return self.q

    def identity(self):
        """
        Returns additive identity, the point at infinity.
        """
        return (0, 0)

    def inverse(self, a):
        """
        Returns inverse (or negative) of a, its reflection across the y axis.
        """
        return (a[0], -a[1])

    def is_element(self, a):
        """
        Tests for equality in the Weierstrauss equation y^2 = x^3 + ax + b
        """
        l = pow(int(a[1]), 2, self.p)
        r = (pow(int(a[0]), 3, self.p) + self.a * a[0] + self.b) % self.p
        return l == r or a == self.identity()

    def random_secret(self):
        """
        Generates random secret in the range [1, self.p - 1]
        """
        return random.randrange(1, self.p)

    def random_element(self, secret=None):
        """
        Generates random element in the field f_q
        """
        if secret is None:
            secret = self.random_secret
        return self.multiply(self.generator(),secret)

    # Mathematical operations

    def add(self, p, q):
        """
        Computes P + Q = R using elliptic curve addition (extending P and Q
        to the point at infinity, then negating it).
        """
        x_p, y_p = p
        x_q, y_q = q

        if (x_p == x_q and y_p == -y_q) or (p == q and y_p == 0):
            return (0,0)
        elif (x_p == x_q and y_p == y_q):
            return self.double(p)
        elif p == (0,0):
            return q
        elif q == (0,0):
            return p

        s = ((y_q - y_p) * inverse(x_q - x_p, self.p)) % self.p
        x_r = (s**2 - x_p - x_q) % self.p
        y_r = (s * (x_p - x_r) - y_p) % self.p
        assert(self.is_element((x_r, y_r)))
        return (x_r, y_r)

    def double(self, p):
        """
        Adding the same point to itself.
        """
        p_x, p_y = p
        if p_x == 0 and p_y == 0:
            return (0,0)
        s = ((3 * p_x**2 + self.a) * inverse(2 * p_y, self.p)) % self.p

        r_x = (s * s - 2 * p_x) % self.p
        r_y = (s * (p_x - r_x) - p_y) % self.p

        return (r_x, r_y)

    def multiply(self, k, n):
        """
        Uses repeated doubling method. n is an integer and k is an elliptic
        curve point. Optimized versions of multiply to be implemented.
        """
        db = bin(n)
        T = self.identity()
        lendb = len(db)
        for i in range(len(db)):
            T = self.double(T)
            if db[i] == str(1):
                T = self.add(T, k)
        assert(self.is_element(T))
        return T

    def bytes(self, a):
        return long_to_bytes(a)

    def element(self, a):
        return bytes_to_long(a)

    def encode(self, data):
        """
        Maps data to the curve using modified Koblitz method.
        """
        if isinstance(data, int):
            data = self.bytes(data)
        tmp_data = bytearray(b'\xff' + data + b'0\xff')
        element = self.element(tmp_data)
        for j in range(1, self.k):
            x_j = (self.k * element + j) % self.p
            y_j = self._solve_for_y(x_j)
            if y_j is not None and self.is_element((x_j, y_j)):
                return ((x_j, y_j))

    def decode(self, pt):
        """
        Returns the actual data encoded in the point. Checks that the
        data is well-formed.
        """
        data = self.bytes(((pt[0] - 1) // self.k) % self.q)
        assert data[0] == 0xff and data[-1] == 0xff
        return data[1:-2]

    def sign(self, secret, data):
        """
        Executes ECDSA for data to be signed.
        """
        s = 0
        e = self._hash(data, self.q - 1)
        while s == 0 or r == 0:
            k = self.random_secret()
            kG_x, _ = self.multiply(self.generator(), k)
            r = kG_x
            s = (inverse(k, self.q) * (e + (secret * r)) % self.q) % self.q
        return (r, s)

    def verify(self, a, data, signature):
        """
        Verifies a's signature on data.
        """
        r, s = signature
        a_x, a_y = a
        assert (1 < r < (self.q - 1)) and (1 < s < (self.q - 1))

        e = self._hash(data, self.q - 1)
        w = inverse(s, self.q)

        u_1 = (e * w) % self.q
        u_2 = (r * w) % self.q
        X_x, X_y = self.add(self.multiply(self.generator(), u_1), \
            self.multiply(a, u_2))

        if X_x == 0 and X_y == 0:
            return False
        return (X_x % self.q) == r

    def _solve_for_y(self, x):
        tmp = (pow(x, 3, self.p) + self.a * x + self.b) % self.p
        # test with jacobi symbol (determines if tmp is a quadratic residue)
        if jacobi(tmp, self.p) == 1:
            return square_root_mod_prime(tmp, self.p)
        else:
            return None

class Test(unittest.TestCase):
    def setUp(self):
        self.group = curvep256()
        self.x0 = PrivateKey(self.group)
        self.x1 = PrivateKey(self.group)
        self.y0 = self.x0.public_key()
        self.y1 = self.x1.public_key()
        self.S, self.T = self.test_operations_basic()

    def test_is_element(self):
        self.assertTrue(self.group.is_element((self.group.g)))
        self.assertTrue(self.group.is_element((self.x1.element)))
        self.assertTrue(self.group.is_element((self.x0.element)))
        for i in range(100):
            self.assertTrue(self.group.is_element((self.group.random_element())))

    def test_operations_basic(self):
        # Values obtained from NIST Routines handbook:
        # http://www.nsa.gov/ia/_files/nist-routines.pdf
        x_s = string_to_long("de2444be bc8d36e6 82edd27e 0f271508 617519b3"
            "221a8fa0 b77cab39 89da97c9")
        y_s = string_to_long("c093ae7f f36e5380 fc01a5aa d1e66659 702de80f"
            "53cec576 b6350b24 3042a256")
        S = (x_s, y_s)

        x_t = string_to_long("55a8b00f 8da1d44e 62f6b3b2 5316212e 39540dc8"
            "61c89575 bb8cf92e 35e0986b")
        y_t = string_to_long("5421c320 9c2d6c70 4835d82a c4c3dd90 f61a8a52"
            "598b9e7a b656e9d8 c8b24316")
        T = (x_t, y_t)
        return (S, T)

    def test_addition(self):
        x_r = string_to_long("72b13dd4 354b6b81 745195e9 8cc5ba69 70349191"
            "ac476bd4 553cf35a 545a067e")
        y_r = string_to_long("8d585cbb 2e1327d7 5241a8a1 22d7620d c33b1331"
            "5aa5c9d4 6d013011 744ac264")
        R = self.group.add(self.S, self.T)
        self.assertEqual(R, (x_r, y_r))

    def test_subtraction(self):
        x_rs = string_to_long("c09ce680 b251bb1d 2aad1dbf 6129deab 837419f8"
            "f1c73ea1 3e7dc64a d6be6021")
        y_rs = string_to_long("1a815bf7 00bd8833 6b2f9bad 4edab172 3414a022"
            "fdf6c3f4 ce30675f b1975ef3")
        R_s = self.group.add(self.S, self.group.inverse(self.T))
        self.assertEqual(R_s, (x_rs, y_rs))

    def test_double(self):
        x_rd = string_to_long("7669e690 1606ee3b a1a8eef1 e0024c33 df6c22f3"
            "b17481b8 2a860ffc db6127b0")
        y_rd = string_to_long("fa878162 187a54f6 c39f6ee0 072f33de 389ef3ee"
            "cd03023d e10ca2c1 db61d0c7")
        R_d = self.group.double(self.S)
        self.assertEqual(R_d, (x_rd, y_rd))

    def test_multiplication(self):
        d = string_to_long("c51e4753 afdec1e6 b6c6a5b9 92f43f8d d0c7a893"
            "3072708b 6522468b 2ffb06fd")
        x_rm = string_to_long("51d08d5f 2d427888 2946d88d 83c97d11 e62becc3"
            "cfc18bed acc89ba3 4eeca03f")
        y_rm = string_to_long("75ee68eb 8bf626aa 5b673ab5 1f6e744e 06f8fcf8"
            "a6c0cf30 35beca95 6a7b41d5")
        R_m = self.group.multiply(self.S, d)
        self.assertEqual(R_m, (x_rm, y_rm))

        # Joint scalar multiply testing
        e = string_to_long("d37f628e ce72a462 f0145cbe fe3f0b35 5ee8332d"
            "37acdd83 a358016a ea029db7")
        x_j = string_to_long("d867b467 92210092 34939221 b8046245 efcf5841"
            "3daacbef f857b858 8341f6b8")
        y_j = string_to_long("f2504055 c03cede1 2d22720d ad69c745 106b6607"
            "ec7e50dd 35d54bd8 0f615275")

        R_j = self.group.add(self.group.multiply(self.S, d), \
            self.group.multiply(self.T, e))
        self.assertEqual(R_j, (x_j, y_j))

    def test_addition_random(self):
        for i in range(50):
            r1 = self.group.random_element()
            r2 = self.group.random_element()
            self.group.add(r1, r2)

    def test_multiplication_random(self):
        for i in range(50):
            k = self.group.random_secret()
            self.group.multiply(self.group.generator(), k)

    def test_encoding(self):
        g = self.group
        msg1 = b"" # 0 character encoding
        msg2 = b"Hello" # multiple character encoding
        msg3 = b"abcdefghijklmnopqrstuvxyzab"
        # Maximum number of chars: 28 + 2 zero padding + leading byte = 31

        self.assertEqual(g.decode(g.encode(msg1)), msg1)
        self.assertEqual(g.decode(g.encode(msg2)), msg2)
        self.assertEqual(g.decode(g.encode(msg3)), msg3)

    def test_ecdsa(self):
        msg = b"Example of ECDSA with P-256"
        g = self.group
        self.assertTrue(g.verify(self.y0.element, msg, \
            g.sign(self.x0.secret, msg)))
        self.assertFalse(g.verify(self.y1.element, msg, \
            g.sign(self.x0.secret, msg)))

    def test_exchange(self):
        self.assertEqual(self.x0.exchange(self.y1), self.y1.exchange(self.x0))
        self.assertEqual(self.x1.exchange(self.y0), self.y0.exchange(self.x1))

if __name__ == "__main__":
    unittest.main()
