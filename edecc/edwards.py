"""
Elliptic curve implementation of the abstract class Group, Curve, and Point.
"""

import math
from elgamal import ElGamal
from group import Group, Point
from ed25519 import ed25519
from modular import ModInt
from utils import b2l, l2b

b = 256 # word size
k = 1000

class EdwardsCurve(Group, ElGamal, object):
    """
    A twisted Edwards curve is described by the equation
    ax^2 + y^2 = c*(1 + dx^2y^2)

    Attributes:
        p(int): Size of the prime field, used in modulo operations.
        r(int): Order of the finite prime field.
        d(int), a(int), c(int): equation parameters

    The EdCurve class generates points, secrets on a twisted Edwards curve.
    """

    def __init__(self, name, p, d, a, r, gx, gy):
        self.name = name
        self.p = p
        self.d = d
        self.a = a
        self.r = r
        self.c = self

        self.one = ModInt(p, 1)
        self.zero = ModInt(p, 0)

        self.base = EdwardsPoint(self, gx, gy)
        self.i = EdwardsPoint(self, self.zero, self.one)

        if not self.base._on_curve():
            raise Exception("Incorrect base point")

    def string(self):
        return self.name

    def point(self):
        return EdwardsPoint(self, ModInt(self.p), ModInt(self.p))

    def secret(self):
        return ModInt(self.p).random_secret()

    def order(self):
        """
        Returns the order of the curve, or the number of unique elements.
        """
        return self.r

    def to_ec_from_tec(self):
        """
        Returns a standard, non-twist Edwards curve from twist curve params.
        """
        newd = ModInt(self.p)
        return EdwardsCurve(self.p, new.div(self.d, self.a), self.a,
                            self.r, self.base.x, self.base.y)

    def to_tec_from_ec(self):
        """
        Returns a twisted Edwards curve from non-twist curve parameters.
        """
        aa = ModInt(self.p)
        return EdwardsCurve(self.p, self.d, aa.mul(self.a, self.a), self.r,
                            self.base.x, self.base.y)

class EdwardsPoint(Point, object):
    def __init__(self, curve, x=ModInt(), y=ModInt()):
        self.c = curve
        self.x = x
        self.y = y

    def string(self):
        return (self.x.v, self.y.v)

    def _on_curve(self):
        """ Determines if this point is a point on its parent curve.

        Args:
            a: point to test

        Returns:
            true if a is a member of the set described by the equation
            cax^2 + y^2 = 1 + dx^2y^2, (twisted Edwards Coordinates.)
            false otherwise.
        """
        x, y = self.x, self.y
        P = self.c.p
        l, r, xx, yy = ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        xx.mul(x, x)
        yy.mul(y, y)
        l.mul(self.c.a, xx).add(l, yy)
        r.mul(self.c.d, xx).mul(r, yy).add(self.c.one, r)
        return l.equal(r)

    def generator(self):
        """"
        Returns base point of the curve.
        """
        self.set(self.c.base)
        return self

    def identity(self):
        """
        Returns additive identity, the neutral element (0, 1)
        """
        self.set(self.c.i)
        return self

    def neg(self, a):
        """
        Returns negation of a, its reflection across the y axis.
        """
        self.x.neg(a.x)
        self.y.set(a.y)
        return self

    def equal(self, p):
        """
        Checks for equality of two points on curve.
        """
        return self.x.equal(p.x) and self.y.equal(p.y)

    def set(self, p):
        """
        Sets point to be equal to p.
        """
        self.x.set(p.x)
        self.y.set(p.y)
        return self

    def add(self, p, q):
        """
        Computes P + Q = R using Edwards elliptic curve addition.
        Complete when a is a square and d is a nonsquare

        x =  (x1y2+x2y1)/(1 + dx1x2y1y2)
        y = (y1y2-ax1x2)/(1 - dx1x2y1y2)
        """
        P = self.c.p
        x1, y1 = p.x, p.y
        x2, y2 = q.x, q.y

        dm = ModInt(P)
        t1, t2 = ModInt(P), ModInt(P)
        nx, dx = ModInt(P), ModInt(P)
        ny, dy = ModInt(P), ModInt(P)

        dm.mul(self.c.d, x1).mul(dm, x2).mul(dm, y1).mul(dm, y2)
        nx.add(t1.mul(x1, y2), t2.mul(x2, y1))
        dx.add(self.c.one, dm)
        ny.sub(t1.mul(y1, y2), t2.mul(x1, x2).mul(self.c.a, t2))
        dy.sub(self.c.one, dm)

        self.x.div(nx, dx)
        self.y.div(ny, dy)

        # x = ((x1.v*y2.v + x2.v * y1.v) * inverse(1 + self.c.d.v * x1.v * x2.v * y1.v * y2.v, P.v))% P.v
        # y =  ((y1.v*y2.v - self.c.a.v * x1.v * x2.v) * inverse(1 - self.c.d.v * x1.v * x2.v * y1.v * y2.v, P.v)) % P.v
        # print("actual", x, y)
        # print("this thing", self.string())
        #assert self._on_curve()
        return self

    def sub(self, p, q):
        """
        Computes p - q by first negating q.
        """
        return self.add(p, q.neg(q))

    def double(self, p):
        """
        Doubles the point p on the curve using its unified addition law.
        """
        return p.add(p, p)

    def multiply(self, P, n):
        """
        Uses repeated doubling method. n is an integer and P is an elliptic
        curve point.
        """
        db = bin(n.v)
        self.set(self.c.i)
        for i in range(len(db)):
            self.double(self)
            if db[i] == str(1):
                self.add(self, P)
        return self

    def multiply_window(self, P, n):
        pass

    def multiply_ladder(self, P, n):
        pass

    def random_element(self, secret=None):
        """
        Generates a random element in the field p.
        """
        if secret is None:
            secret = self.c.secret()
        G = self.c.point()
        G.set(self.generator())
        self.multiply(G, secret)
        #assert self._on_curve()
        return self

    def _special_pt(self, a):
        """
        For inverted coordinate system; checks if a is one of four special
        points.
        """

    def to_ep(self, a):
        return a

    def from_ep(self, a):
        return a

    def encode(self, data):
        """
        Maps data to the curve using modified Koblitz method.
        here, k is the number of tries.

        Generates a point on the Edwards curve.
        """
        P = self.c.p
        K = ModInt(P, k)
        if isinstance(data, int):
            data = l2b(data)
        tmp_data = bytearray(b'\xff' + data + b'0\xff')
        if len(data) == 0:
            tmp_data = bytearray(b'\xff' + b'\xff')
        element = ModInt(P, b2l(tmp_data))
        for j in range(k):
            self.y.mul(K, element).add(self.y, ModInt(P, j))
            pt = self._solve_for_x(self.y)
            if pt is not None:
                return pt

    def _solve_for_x(self, y):
        """
        Given the y-coordinate, return the x coordinate of the corresponding
        point on the curve. Because this is only used for encoding plaintext
        messages, either root is okay (as long as it is a point on the curve).

        The Jacobi symbol is used to quickly determine if xx is a quadratic
        residue.

        x^2 = (y^2 - 1) / (dy^2 + 1)
        """
        P = self.c.p
        xx, yy, num, denom, x = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        yy.mul(y, y)
        num.sub(yy, self.c.one)
        denom.add(denom.mul(yy, self.c.d), self.c.one)
        xx.div(num, denom)
        if self.x.jacobi(xx) == 1:
            x.sqrt(xx)
            newpt = EdwardsPoint(self.c, x, y)
            if not newpt._on_curve():
                x.sub(P, self.x)
            return self.from_ep(newpt)

    def decode(self, pt):
        """
        Returns the actual data encoded in the point. Checks that the
        data is well-formed.

        Issues: zero padding doesn't work sometimes?
        """
        data = l2b((pt.to_ep(pt).y.v - 1) // k)
        return data[1:-2]

Group.register(EdwardsCurve)
Point.register(EdwardsPoint)

Group.__subclasscheck__(EdwardsCurve)
Point.__subclasscheck__(EdwardsPoint)
Group.__instancecheck__(EdwardsCurve)
Point.__instancecheck__(EdwardsPoint)