import edwards
from modular import ModInt
from group import Group, Point

#TODO
#   - test:
#       * clear_denom
#       * other doubling algorithm
#   - edge/special points?

class projEdwardsCurve(edwards.EdwardsCurve, object):
    """
    The ProjectiveCurve class instantiates curves in Projective coordinates.
    An Edwards curve in Projective coordinates describes the set of points in the
    equation (x^2 + ay^2) * z^2 = z^4 + dx^2y^2.

    This class executes arithmetic operations in Projective coordinates,
    and inherits the rest of the curve functions from the EdwardsCurve class.

    To convert from standard Edwards coordinates to projective coordinates,
    compute (X, Y, 1).

    The addition formulas are still strongly unified (so they can
    be used to double a point).

    Attributes:
        p (int): Order of the finite prime field that the curve is defined over.
        a, d (int): Parameters of the equation.

    """

    def __init__(self, ed):
        self.name = "Edwards projective"
        self.c = ed
        self.a = ed.a
        self.d = ed.d
        self.p = ed.p
        self.r = ed.r

        self.zero = ed.zero
        self.one = ed.one

        self.base = self.point().from_ep(ed.base)
        self.i = projEdwardsPoint(self, self.zero, self.one, self.one)
        if not self.base._on_curve():
            raise Exception("Incorrect base point")

    def point(self):
        """
        Generates a new projective Edwards point.
        """
        return projEdwardsPoint(self, ModInt(self.p), ModInt(self.p), ModInt(self.p))

class projEdwardsPoint(edwards.EdwardsPoint, object):
    def __init__(self, curve, x=ModInt(), y=ModInt(), z=ModInt()):
        self.c = curve
        self.x = x
        self.y = y
        self.z = z

    def string(self):
        return (self.x.v, self.y.v, self.z.v)

    def _on_curve(self):
        """
        Returns if a is a member of the set described by the equation
        (ax^2 + y^2) * z^2 = z^4 + dx^2y^2, that is, in
        projective twisted Edwards Coordinates.
        """
        x, y, z = self.x, self.y, self.z
        P = self.c.p
        l, r, xx, yy, zz = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        m = ModInt(P)

        xx.mul(x, x)
        yy.mul(y, y)
        zz.mul(z, z)
        l.mul(self.c.a, xx).add(l, yy).mul(l, zz)
        m.mul(self.c.d, xx).mul(m, yy)
        r.exp(zz, ModInt(P, 2)).add(r, m)
        return l.equal(r)

    def identity(self):
        """
        The identity element is (0, 1, 1).
        """
        self.set(self.c.i)
        return self

    def neg(self, a):
        """
        The negative of element (x, y, z) is (-x, y, z)
        """
        self.x.neg(a.x)
        self.y.set(a.y)
        self.z.set(a.z)
        return self

    def equal(self, p):
        return self.x.equal(p.x) and self.y.equal(p.y) and self.z.equal(p.z)

    def set(self, p):
        self.x.set(p.x)
        self.y.set(p.y)
        self.z.set(p.z)
        return self

    def add(self, p, q):
        """
        Adds p and q, both expressed as projective coordinates.
        Computational cost: 10M + 1S + 2D + 7add.
        """
        x1, y1, z1 = p.x, p.y, p.z
        x2, y2, z2 = q.x, q.y, q.z
        P = self.c.p
        A, B, C, D, E, F, G  = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        t1, t2, t3, t4 = ModInt(P), ModInt(P), ModInt(P), ModInt(P)

        # print("\nAdding: ", "P", p.string(), "q", q.string())
        if z1.equal(self.c.zero) or z2.equal(self.c.zero):
            raise Exception("Point", p.string(), "or point", q.string(), "not representable")

        # elif z1.equal(self.c.one) and z2.equal(self.c.one):
        #     A.set(self.c.one)
        # elif z1.equal(self.c.one):
        #     A.set(z2)
        # elif z2.equal(self.c.one):
        #     A.set(z1)
        # else:
        A.mul(z1, z2)
        B.mul(A, A)
        C.mul(x1, x2)
        D.mul(y1, y2)
        E.mul(self.c.d, C).mul(E, D)
        F.sub(B, E)
        G.add(B, E)

        t1.add(x1, y1)
        t2.add(x2, y2)
        t4.mul(t1, t2).sub(t4, C).sub(t4, D)
        self.x.mul(A, F).mul(self.x, t4)
        t3.sub(D, t3.mul(self.c.a, C))
        self.y.mul(A, G).mul(self.y, t3)
        self.z.mul(F, G)
        #assert self._on_curve()
        return self

    def double(self, p):
        """
        Doubles p, expressed as projected coordinate.
        Computational cost: 3M + 4S + 1d + 7add
        """
        x, y, z = p.x, p.y, p.z
        # print("\nDoubling", p.string())
        P = self.c.p
        two = ModInt(P, 2)
        B, C, D, E, F, H, J = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)

        B.exp(B.add(x, y), two)
        C.mul(x, x)
        D.mul(y, y)
        E.mul(self.c.a, C)
        F.add(E, D)
        if z.equal(self.c.one):
            H.set(self.c.one)
        else:
            H.mul(z, z)
        J.sub(F, J.mul(two, H))
        self.x.sub(B, C).sub(self.x, D).mul(self.x, J)
        self.y.mul(F, self.y.sub(E, D))
        self.z.mul(F, J)
        #assert self._on_curve()
        return self

    def clear_denom(self, p, q):
        """
        Alternate addition approach with projective twisted edwards coordinates.
        Use substitution: x = sqrt(a) * xbar, y = ybar
        Computational cost: 10M + 1S + 3D + 7add
        """
        x_1, y_1, z_1 = p
        x_2, y_2, z_2 = q

        A = z_1 * z_2
        B = self.a * A * A
        H = self.a * A
        C = x_1 * x_2
        D = y_1 * y_2
        E = self.d * C * D
        F = B - E
        G = B + E

        x_3 = (H * F * ((x_1 + y_1) * (x_2 + y_2) - C - D)) % self.p
        y_3 = (H * G * (D - C)) % self.p
        z_3 = (F * G) % self.p

        return (x_3, y_3, z_3)

    def from_ep(self, a):
        """
        Converts a point in Edwards coordinates to an equivalent point in
        projective coordinates, with respect to the parent curve.
        """
        assert isinstance(a, edwards.EdwardsPoint)
        self.x.set(a.x)
        self.y.set(a.y)
        self.z.set(a.c.one)
        return self

    def to_ep(self, a):
        """
        Converts a, a point in projective Edwards coordinates, to a point in
        standard Edwards coordinates using the transformation
        x = X/Z, y = Y/Z
        """
        X, Y, Z = a.x, a.y, a.z
        x, y = ModInt(a.c.p), ModInt(a.c.p)
        x.div(a.x, a.z)
        y.div(a.y, a.z)
        ed = self.c.c.point()
        ed.x.set(x)
        ed.y.set(y)
        return ed

Group.register(projEdwardsCurve)
Point.register(projEdwardsPoint)

Group.__subclasscheck__(projEdwardsCurve)
Point.__subclasscheck__(projEdwardsPoint)
Group.__instancecheck__(projEdwardsCurve)
Point.__instancecheck__(projEdwardsPoint)