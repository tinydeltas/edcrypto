import edwards
from modular import ModInt
from group import Point, Group

#TODO
#   - test:
#       - fast addition/doubling
#       - detect places where z = 1

class extEdwardsCurve(edwards.EdwardsCurve, object):
    """
    The ExtendedCurve class instantiates curves in Extended coordinates.
    An Edwards curve in Extended coordinates describes the set of points in the
    equation (ax^2 + y^2) = 1 + dx^2y^2.
    However, the Extended representation includes several more parameters:
    (X, Y, T, Z), T = XY/Z
    This corresponds to the extended affine point (X/Z, Y/Z, T/Z).
    However, this system of extended coordinates only apply when a = -1,
    as in the case of ed25519.
    A separate set of extended coordinates applies in other cases.

    This class performs arithmetic operations (e.g. addition)
    in extended coordinates,
    and inherits the rest of the curve functions from the EdwardsCurve class.

    Extended coordinates preserve the strong unification property of
    Edwards coordinates.

    To convert from Twisted Edwards coordinates to Extended coordinates,
    use the maps
    E(x, y) --> P(X, Y, Z) = (x, y, t) == (x, y, xy)
    P(X, Y, Z) --> Ex(X, Y, Z, T) =  (x, y, xy, 1)

    To convert from the Edwards curve to a more compact form
    (ideal for addition in Extended coordinates):
    ax^2 + y2 = 1 + dx^2y^2 --> -x^2+y^2 = 1 + (-d/a)x^2y^2
    if -a is a square in K.

    Attributes:
        p (int): Order of the finite prime field that the curve is defined over.
        a, d (int): Parameters of the equation.

    Source:
    Twisted Edwards Curves revisited. http://eprint.iacr.org/2008/522
    """

    def __init__(self, ed):
        self.name = "Edwards extended"
        self.c = ed
        self.a = ed.a
        self.d = ed.d
        self.p = ed.p
        self.r = ed.r

        self.zero = ed.zero
        self.one = ed.one

        self.base = self.point().from_ep(ed.base)
        self.i = extEdwardsPoint(self, self.zero, self.one, self.zero, self.one)
        if not self.base._on_curve():
            raise Exception("Incorrect base point")

    def point(self):
        """
        Generates a new extended Edwards point.
        """
        return extEdwardsPoint(self, ModInt(self.p), ModInt(self.p),
                                ModInt(self.p), ModInt(self.p))

class extEdwardsPoint(edwards.EdwardsPoint, object):
    def __init__(self, curve, x=ModInt(), y=ModInt(), t=ModInt(), z=ModInt()):
        self.c = curve
        self.x = x
        self.y = y
        self.t = t
        self.z = z

    def string(self):
        return (self.x.v, self.y.v, self.t.v, self.z.v)

    def _on_curve(self):
        """
        Determines if the Extended point is on the curve by first
        converting it into an Edwards point, and then checking the usual way.
        """
        ed = self.to_ep(self)
        x, y = ed.x, ed.y
        p = self.c.p
        xx, yy, l, r = ModInt(p), ModInt(p), ModInt(p), ModInt(p)
        xx.mul(x, x)
        yy.mul(y, y)
        l.mul(self.c.a, xx).add(l, yy)
        r.mul(self.c.d, xx).mul(r, yy).add(self.c.one, r)
        return l.equal(r)

    def identity(self):
        """
        The identity element is (0, 1, 0, 1)
        """
        self.set(self.c.i)
        return self

    def neg(self, a):
        """
        In extended coordinates, the negative of an element
        (x, y, t, z) is (-x, y, -t, z).
        """
        self.x.neg(a.x)
        self.y.set(a.y)
        self.t.neg(a.t)
        self.z.set(a.z)
        return self

    def equal(self, p):
         return self.x.equal(p.x) and self.y.equal(p.y) and self.z.equal(p.z) and self.t.equal(p.t)

    def set(self, p):
        self.x.set(p.x)
        self.y.set(p.y)
        self.t.set(p.t)
        self.z.set(p.z)
        return self

    def add(self, p, q):
        """
        Unified addition in Extended coordinates.
        Adds p and q, both expressed as projective coordinates.
        Computational cost:  9M + 1*a + 1*d + 7add
        """
        x1, y1, t1, z1 = p.x, p.y, p.t, p.z
        x2, y2, t2, z2 = q.x, q.y, q.t, q.z
        P = self.c.p
        A, B, C, D, E, F, G, H = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        T1, T2 = ModInt(P), ModInt(P)

        A.mul(x1, x2)
        B.mul(y1, y2)
        C.mul(t1, t2).mul(C, self.c.d)
        D.mul(z1, z2)

        T1.add(x1, y1)
        T2.add(x2, y2)
        E.mul(T1, T2).sub(E, A).sub(E, B)
        F.sub(D, C)
        G.add(D, C)
        H.sub(B, H.mul(self.c.a, A))

        self.x.mul(E, F)
        self.y.mul(G, H)
        self.t.mul(E, H)
        self.z.mul(F, G)
        #assert self._on_curve()
        return self

    def double(self, p):
        """
        Dedicated doubling in Extended coordinates.
        Independent of curve constant D .
        Doubles p, expressed as projected coordinate.
        Computational cost:  4M + 4S + 1*a + 6add + 1*2.
        """
        x, y, t, z = p.x, p.y, p.t, p.z
        P, two = self.c.p, ModInt(self.c.p, 2)
        assert not z.equal(self.c.zero)
        A, B, C, D, E, F, G, H = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)

        A.mul(x, x)
        B.mul(y, y)
        C.mul(two, C.mul(z, z))
        D.mul(self.c.a, A)
        E.mul(E.add(x, y), E.add(x, y))
        E.sub(E, A)
        E.sub(E, B)

        G.add(D, B)
        H.sub(D, B)
        F.sub(G, C)

        self.x.mul(E, F)
        self.y.mul(G, H)
        self.t.mul(E, H)
        self.z.mul(F, G)
        #assert self._on_curve()
        return self

    def double_fast(self, p):
        """
        Doubling when Z = 1
        Cost: 3M + 4S + 1*a + 7add + 1*2
        """
        x, y, z = p.x, p.y, p.z
        P , two = self.c.p, ModInt(self.c.p, 2)
        A, B, D, E, G, H = ModInt(P),  ModInt(P),  ModInt(P),  ModInt(P),  ModInt(P),  ModInt(P)
        t1, t2, t3 = ModInt(P),  ModInt(P),  ModInt(P)

        A.mul(x, x)
        B.mul(y, y)
        D.mul(self.c.a, A)
        t1.add(x1, y1)
        E.mul(t1, t1).sub(E, A).sub(E, B)
        G.add(D, B)
        H.sub(D, B)

        t2.sub(G, two)
        self.x.mul(E, t2)
        self.y.mul(G, H)
        self.t.mul(E, H)
        t3.mul(two, G)
        self.z.mul(G, G).sub(self.z, t3)

        return self

    def add_fast(self, p, q):
        """
        Fast addition in Extended coordinates is possible when the curve is
        transformed into an equivalent Twisted Edwards curve by the map
        (x, y) -> (x/sqrt(-a), y)
        Thus addition is performed on the curve
        (-x^2 + y^2) = 1 + d'x^2y^2 for d' = -d/a
        Computational cost: 8M + 1D (7M + 1D when z2 = 1)
        """
        x1, y1, t1, z1 = p.x, p.y, p.t, p.z
        x2, y2, t2, z2 = q.x, q.y, q.t, q.z

        xi, yi, zi = self.add_fast_test(p, q)

        P, two = self.c.p, ModInt(self.c.p, 2)

        A, B, C, D, E, F, G, H = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)

        t1, t2, t3, t4 = ModInt(P), ModInt(P), ModInt(P), ModInt(P)

        t1.sub(y1, x1)
        t2.sub(y2, x2)
        t3.add(y1, x2)
        t4.add(y2, x2)

        A.mul(t1, t2)
        B.mul(t3, t4)
        C.mul(self.c.d, two).mul(C, t1).mul(C, t2)
        D.mul(two, z1).mul(D, z2)
        E.sub(B, A)
        F.sub(D, C)
        G.add(D, C)
        H.add(B, A)

        self.x.mul(E, F)
        self.y.mul(G, H)
        self.t.mul(E, H)
        self.z.mul(F, G)

        assert self._on_curve()
        return self

    def add_fast_test(self, p, q):

        x1, y1, t1, z1 = p.x.v, p.y.v, p.t.v, p.z.v
        x2, y2, t2, z2 = q.x.v, q.y.v, q.t.v, q.z.v

        A = ((y1 - x1) * (y2 - x2))
        B = ((y1 + x2) * (y2 + x2))
        C = self.d * 2 * t1 * t2
        D = 2 * z1 * z2
        E = B - A
        F = D - C
        G = D + C
        H = B + A

        x_r = (E * F) % self.p
        y_r = (G * H) % self.p
        t_r = (E * H) % self.p
        z_r = (F * G) % self.p

        p = (x_r, y_r, t_r, z_r)
        return p

    def to_ep(self, a):
        """
        Converts a, a point in extended Edwards coordinates, to a point in
        standard Edwards coordinates.
        Alternative method: ext -> projective -> regular
        """
        X, Y, T, Z = a.x, a.y, a.t, a.z
        x, y = ModInt(a.c.p), ModInt(a.c.p)
        t1, t2 = ModInt(a.c.p), ModInt(a.c.p)

        x.div(X, Z)
        y.div(Y, Z)

        ed = self.c.c.point()
        ed.x.set(x)
        ed.y.set(y)
        assert t1.mul(x, y).equal(t2.div(T, Z))
        return ed

    def from_ep(self, a):
        """
        Converts an Edwards point to the equivalent point on the
        extended coordinate system.
        """
        assert isinstance(a, edwards.EdwardsPoint)
        z = ModInt(a.c.p)
        self.x.set(a.x)
        self.y.set(a.y)
        self.t.set(z.mul(a.x, a.y))
        self.z.set(a.c.one)
        return self

Group.register(extEdwardsCurve)
Point.register(extEdwardsPoint)

Group.__subclasscheck__(extEdwardsCurve)
Point.__subclasscheck__(extEdwardsPoint)
Group.__instancecheck__(extEdwardsCurve)
Point.__instancecheck__(extEdwardsPoint)

