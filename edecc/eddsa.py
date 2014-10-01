"""
Performs ed25519, a variant of the ECDSA algorithm specifically
for the Edwards curve ed25519, on Edwards public/private key pairs.
"""

from edwards import EdwardsCurve, EdwardsPoint
from utils import *
from modular import ModInt

b = 256

class edwardsPublicKey(EdwardsCurve):
    def __init__(self, group, element):
        self.group = group
        self.element = element

    def _hash(self, msg):
        """
        Custom hash function for edDSA: public keys are 2b bits long.
        Uses a cryptographically secure sha512 function.
        """
        hasher = SHA512.new()
        hasher.update(bytes(msg, encoding='utf-8'))
        return hasher.digest()

    def sum_bits(self, bits):
        """
        Faster version of sum_bits_a
        """
        return ''.join([chr(sum([bits[i * 8 + j] << j \
                for j in range(8)])) for i in range(b//8)])

    def bit(self, hashed, i):
        """
        Given a (hashed) string h, return the value of the bit at index i.
        """
        hashed = str(hashed)
        return (ord(hashed[i//8]) >> (i % 8)) & 1

    def _Hint(self, m):
        """
        Given a message m, hash it using the SHA512 method and convert that
        value into its (2b) bit-array representation.
        """
        h = self._hash(m)
        return sum(pow(2, i, self.group.p.v) * self.bit(h, i) for i in range(2*b))

    def _decodeint(self, s):
        """
        Assume int s has been encoded by the encodeint method above.
        Converts from a 256-bit-array to a base-10 integer.
        """
        return sum(pow(2, i, self.group.p.v) * self.bit(s, i) for i in range(0, b))

    def _decodepoint(self, s):
        """
        Given the bit-array representation s (actually the encoded y-coordinate
        of some point), convert the first (b-1) bits into a base-10 integer,
        recover x, and change the sign if necessary.
        """
        p = self.group.p
        d = self.group.d
        yy, xx, x, num, denom, test = ModInt(p), ModInt(p), ModInt(p), ModInt(p), ModInt(p), ModInt(p)

        y = ModInt(p, sum(pow(2, i, p.v) * self.bit(s, i) for i in range(0, b-1)))
        yy.mul(y, y)
        num.sub(yy, self.group.one)
        denom.add(denom.mul(yy, self.group.d), self.group.one)
        xx.div(num, denom)
        x.exp(xx, ModInt(p, ((p.v + 3)//8)))

        if not (test.mul(x, x).sub(test, xx)).equal(self.group.zero):
            I = ModInt(p)
            I.exp(ModInt(p, 2), ModInt(p, (p.v-1)//4))
            x.mul(x, I)
        if x.v % 2 != 0:
            x.sub(p, x)
        if x.v & 1 != self.bit(s, b - 1):
            x.sub(p, x)

        P = self.group.point()
        ed = self.group.c.point()
        ed.x.set(x)
        ed.y.set(y)
        if not ed._on_curve():
            raise Exception("decoding point that is not on curve")
        return P.from_ep(ed)

    def verify(self, m, s):
        """
        Verifies a message m using the edDSA algorithm.
        Checks that R, A, and s values have been recreated correctly by
        making another commitment (hashhing r, pk, and m together), as well as
        checking that
        sB = rB + H(Rpub, Apub, M)aB = R + H(Rpub, Apub, M)A.
        """
        r, s = s
        pk = self.element
        p = self.group.p
        if (len(r) + len(s)) != b//4:
            raise Exception("signature length is wrong")
        if len(pk) != b//8:
            raise Exception("public-key length is wrong")
        R = self._decodepoint(r)
        A = self._decodepoint(pk)
        S = self._decodeint(s)
        h = self._Hint("%s%s%s" % (r, pk, m))
        # print("R", R.string())
        # print("A", A.string())
        # print("S", S)
        # print("h", h)
        elem1 = self.group.point()
        return elem1.random_element(ModInt(p, 8*S)) == elem1.add(elem1.multiply(R, ModInt(p, 8)), elem1.multiply(A, ModInt(p, 8*h)))

class edwardsPrivateKey(edwardsPublicKey):
    def __init__(self, group, secret = None):
        """
        Generates a random private/public key pair in the field f_q
        Takes a secret (a random integer in the field)
        """
        if secret is None:
            secret = group.secret().v
        h = str(self._hash(self._encodeint(secret)))
        p = group.p.v
        a_1 = pow(2, b-2, p)
        a_2 = sum(pow(2, i, p) * self.bit(h, i) for i in range(3, b-2))
        a = ModInt(p, a_1 + a_2)
        P = group.point()
        element = self._encodepoint(P.to_ep(P.random_element(a)))
        edwardsPublicKey.__init__(self, group, element)
        self.secret = secret

    def _encodeint(self, y):
        """
        Given an integer, returns a string of ascii characters representing y.
        y is first converted to its (little-endian) bit-array representation.
        (e.g. 4 --> [0, 0, 1, 0.....0] for array of size 256), which is
        then converted to an ascii-string-literal representation.
        """
        bits = [(y >> i) & 1 for i in range(b)]
        return self.sum_bits(bits)

    def _encodepoint(self, P):
        """
        Given a point on the curve, returns a string of ascii characters that
        represent the point.
        First convert the y coordinate to its b - 1 bit-array representation,
        with the last bit being the sign of the x coordinate.
        It is then converted to an ascii-string-literal representation.
        """
        x, y = P.x.v, P.y.v
        bits = [(y >> i) & 1 for i in range(b - 1)] + [x & 1]
        return self.sum_bits(bits)

    def sign(self, m):
        """
        Signs a message m using the edDSA algorithm. Described here:
        http://ed25519.cr.yp.to/ed25519-20110926.pdf

        Args:
            - (array) sk, the encoded version of the (int) private key
            - (array) pk, the encoded version of the (int) public key
            - (string) m, the message to be signed

        Returns:
            - (array) r, the encoded version of the point multiplied by the last
            (b, 2b-1) bits of the bit-array determined by the hash of the
            secret key.
            - (array) s, the encoded version of the int determined by
            calculating the commitment (hashing)
        """
        p = self.group.p.v
        pk = self.element
        h = str(self._hash(self._encodeint(self.secret)))
        a_1 = pow(2, b-2, p)
        a_2 = sum(pow(2, i, p) * self.bit(h, i) for i in range(3, b-2))
        a = a_1 + a_2
        r = self._Hint(''.join([h[i] for i in range(b//8, b//4)]) + m)
        R = self.group.point()
        R = R.to_ep(R.random_element(ModInt(p, r)))
        S = (r + self._Hint("%s%s%s" % (self._encodepoint(R), pk, m)) * a) % self.group.r.v
        # print("a", a)
        # print("r", r)
        # print("R", R.string(), R._on_curve(), type(R), R.c.name)
        # print("S",S)
        return (self._encodepoint(R), self._encodeint(S))

    def public_key(self):
        return edwardsPublicKey(self.group, self.element)


