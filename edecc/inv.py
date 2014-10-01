"""
InvertedCurve module implements the Inverted representation of a
twisted Edwards interface and its associated unit tests.
"""

import edwards
from modular import ModInt
from group import Group, Point

# todo
#   - test:
#       * tripling 1, 2
#       * other doubling algorithm
#       * clear denom 1, 2

class invEdwardsCurve(edwards.EdwardsCurve, object):
    """
    The InvertedCurve class instantiates curves in Inverted coordinates.
    An Edwards curve in inverted coordinates describes the set of points in the
    equation (x^2 + ay^2) * z^2 = dz^4 + x^2y^2.
    This class executes arithmetic operations in Inverted coordinates,
    and inherits the rest of the curve functions from the EdwardsCurve class.

    To convert from standard Edwards coordinates to inverted coordinates,
    compute (YZ, XZ, XY) (with Edwards coordinates in (X, Y, 1) form).

    This representation does not cover the points (0,+-1) and (+-1,0)
    since xyz != 0 .
    In standard Edwards coordinate, these correspond to points of infinity at
    (0, 1, 1), (0, -1, 1), (1, 0, 1), and (-1, 0, 1)
    In inverted coordinates, the points at infinity are
    (1, 0, 0), (-1, 0, 0), (0, 1, 0), and (0, -1, 0)
    Thus, we must check for z = 0 before adding to the point to another, and
    check for xy = 0 before converting from standard edwards coordinates
    (As you may lose distinction between points.)

    However, the addition formulas are still strongly unified (so they can
    be used to double a point).

    Attributes:
        p (int): Order of the finite prime field that the curve is defined over.
        a, d (int): Parameters of the equation.
    """

    def __init__(self, ed):
        self.name = "Edwards inverted"
        self.c = ed
        self.a = ed.a
        self.d = ed.d
        self.p = ed.p
        self.r = ed.r

        self.zero = ed.zero
        self.one = ed.one
        self.negone = ModInt(self.p, -1)

        self.s1 = invEdwardsPoint(self.c, self.one, self.zero, self.zero)
        self.s2 = invEdwardsPoint(self.c, self.negone, self.zero, self.zero)
        self.s3 = invEdwardsPoint(self.c, self.zero, self.negone, self.zero)
        self.s4 = invEdwardsPoint(self.c, self.zero, self.one, self.zero)

        self.base = self.point().from_ep(ed.base)
        self.i = invEdwardsPoint(self, self.one, self.zero, self.zero)
        if not self.base._on_curve():
            raise Exception("Incorrect base point")

    def point(self):
        """
        Generates a new inverted Edwards point.
        """
        return invEdwardsPoint(self, ModInt(self.p), ModInt(self.p), ModInt(self.p))

class invEdwardsPoint(edwards.EdwardsPoint, object):
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
        (x^2 + ay^2) * z^2 = dz^4 + x^2y^2, that is, in
        inverted twisted Edwards Coordinates.
        """
        x, y, z = self.x, self.y, self.z
        P = self.c.p
        l, r, xx, yy, zz = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        m = ModInt(P)

        xx.mul(x, x)
        yy.mul(y, y)
        zz.mul(z, z)

        l.mul(self.c.a, yy).add(l, xx).mul(l, zz)
        m.mul(xx, yy)
        r.exp(zz, ModInt(P, 2)).mul(r, self.c.d).add(r, m)
        return l.equal(r)

    def identity(self):
        """
        Identity after converting from (0, 1) => (0, 1, 1) => (1, 0, 0)
        """
        self.set(self.c.i)
        return self

    def neg(self, a):
        """
        Corresponding inverse.
        -(x, y, z) = (-x, y, z)
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

    def _special_pt(self, a):
        """
        Checks if a is one of four special points on curve:
        (1, 0, 0), (-1, 0, 0), (0, 1, 0) or (0, -1, 0).
        """
        if a.equal(self.c.s1) or a.equal(self.c.s2) or a.equal(self.c.s3) or a.equal(self.c.s4):
            return True
        return False

    def add(self, p, q):
        """
        Adds p and q, both expressed as inverted coordinates.
        Computational cost: 9M + 1S + 2D + 7add
        Also supports mixed addition: when z_1 is 1, A is simply z_1.

        Note that this also supports the addition of special points:
        1) If z_1 or z_2 = 0,
        p + q = (x_1 * x_2 - y_1 * y_2, x_2 * y_1 + x_1 * y_2, z_1 + z_2)
        2) If I = 0, and y_2 * z_1 = y_1 * z_2, the sum is (1, 0, 0).
        3) If I = 0, and y_2 * z_1 = -y_1 * z_2, the sum is (-1, 0, 0).
        4) If H = 0, and y_2 * z_1 = -x_1 * z_2, the sum is (0, 1, 0).
        5) If H = 0, and y_2 * z_1 = x_1 * z_2, the sum is (0, -1, 0).
        """
        x1, y1, z1 = p.x, p.y, p.z
        x2, y2, z2 = q.x, q.y, q.z
        P, zero = self.c.p, self.c.zero
        A, B, C, D, E, H, I = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        t1, t2, t3, t4, t5, t6 = ModInt(P), ModInt(P), ModInt(P), ModInt(P),  ModInt(P), ModInt(P)

        if z1.equal(self.c.one) and z2.equal(self.c.one):
            A.set(self.c.one)
        elif z1.equal(self.c.one):
            A.set(z2)
        elif z2.equal(self.c.one):
            A.set(z1)
        else:
            A.mul(z1, z2)

        B.mul(self.c.d, A).mul(B, A)
        C.mul(x1, x2)
        D.mul(y1, y2)
        E.mul(C, D)
        H.sub(C, H.mul(self.c.a, D))
        t1.add(x1, y1)
        t2.add(x2, y2)
        I.mul(t1, t2).sub(I, C).sub(I, D)
        t3.mul(y2, z1)
        if H.equal(zero) and t3.equal(t4.mul(t4.neg(x1), z2)):
            self.set(self.c.s4)
        elif H.equal(zero) and t3.equal(t4.mul(x1, z2)):
            self.set(self.c.s3)
        elif I.equal(zero) and t3.equal(t4.mul(y1, z2)):
            self.set(self.c.s1)
        elif I.equal(zero) and t3.equal(t4.mul(t4.neg(y1), z2)):
            self.set(self.c.s2)
        elif z1.equal(zero) or z2.equal(zero) and self._special_pt(p) or self._special_pt(q):
            t5.mul(x2, y1)
            t6.mul(x1, y2)
            self.x.sub(C, D)
            self.y.add(t5, t6)
            self.z.add(z1, z2)
        else:
            self.x.add(E, B).mul(self.x, H)
            self.y.sub(E, B).mul(self.y, I)
            self.z.mul(A, H).mul(self.z, I)
        #assert self._on_curve()
        return self

    def double(self, p):
        """
        Doubling in inverted coordinates.
        Computational cost: 3M + 3S + 1*a + 6add.
        """
        x, y, z = p.x, p.y, p.z
        P, zero, two = self.c.p, self.c.zero, ModInt(self.c.p, 2)
        A, B, U, C, D, E = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        t1, t2, t3, t4, xyz = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)

        A.mul(x, x)
        B.mul(y, y)
        U.mul(self.c.a, B)
        C.add(A, U)
        D.sub(A, U)
        t1.add(x, y)
        E.mul(t1, t1).sub(E, A).sub(E, B)
        xyz.mul(x, y).mul(xyz, z)
        t3.mul(y, z)
        if xyz.equal(zero) and t3.equal(t3):
            self.set(self.c.s1)
        elif xyz.equal(zero) and t3.equal(t4.mul(t4.neg(y), z)):
            self.set(self.c.s2)
        elif xyz.equal(zero) and t3.equal(t4.mul(t4.neg(x), z)):
            self.set(self.c.s4)
        elif xyz.equal(zero) and t3.equal(t4.mul(x, z)):
            self.set(self.c.s3)
        elif z.equal(zero):
            self.x.sub(A, B)
            self.y.mul(x, y).mul(self.y, two)
            self.z.mul(two, z)
        else:
            self.x.mul(C, D)
            t2.mul(z, z).mul(t2, self.c.d).mul(t2, two)
            self.y.sub(C, t2).mul(self.y, E)
            self.z.mul(D, E)
        #assert self._on_curve()
        return self

    def triple(self, p):
        """
        Triples p (inverted coordinate).
        This formula uses 9M + 4S + 1D + 10a.
        """
        x, y, z = p

        A = x * x
        B = y * y
        C = z**2
        D = A + B
        E = 4 * (D - self.d * C)
        H = 2*D * (B - A)
        P = D**2 - A * E
        Q = D**2 - B * E

        x_r = (H + Q) * Q * x
        y_r = (H - P) * P * y
        z_r = P * Q * z

        p = (x_r, y_r, z_r)
        assert self.is_element(p)
        return p

    def triple_i(self, p):
        """
        Triples p (an inverted Edwards coordinate).
        This formula uses 7M + 7S + 1D + 17a, and is faster when
        S/M is small.
        If p is a special point, we have 3p = (x, -y, 0)
        """
        x, y, z = p

        if x*y*z == 0:
            return (x, -y, 0)

        A = x**2
        B = y**2
        C = z**2
        D = A + B
        E = 4 * (D - self.d * C)
        H = 2* D * (B - A)
        P = D**2 - A * E
        Q = D**2 - B * E

        x_r = (H + Q) * ((Q + x)**2 - Q**2 - A)
        y_r = 2 * (H - P) * P * y
        z_r = P * ((Q + z)**2 - Q**2 - C)

        p = (x_r, y_r, z_r)
        #assert self.is_element(p)
        return p

    def clear_denom_i1(self, p, q):
        """
        Alternate addition approach for inverted coordinates.
        Computational cost: 9M + 1s + 3d + 7add.
        """
        x_1, y_1, z_1 = p
        x_2, y_2, z_2 = q

        A = z_1 * z_2
        B = self.d * A * A
        C = x_1 * x_2
        D = y_1 * y_2
        E = self.a * C * D
        H = C - D
        I = (x_1 + y_1) * (x_2 + y_2) - C - D

        x_3 = (E + B) * H
        y_3 = (E - B) * I
        z_3 = self.a * A * H * I

        return (x_3, y_3, z_3)

    def clear_denom_i2(self, p):
        """
        Alternate method for addition in inverted coordinates.
        Computational cost: 3M + 4S + 3D + 5add
        """
        x_1, y_1, z_1 = p

        A = x_1 * x_1
        B = y_1 * y_1
        C = A + B
        D = A - B
        E = pow((x_1 + y_1), 2) - C
        F = self.a * C

        x_r = F * D
        y_r = E * (F - 2 * self.d * z_1 * z_1)
        z_r = self.a * D * E

        return (x_r, y_r, z_r)

    def from_ep(self, a):
        """
        Converts a point in Edwards coordinates to an equivalent point in
        inverted coordinates, with respect to the parent curve.
        """
        assert isinstance(a, edwards.EdwardsPoint)
        one = self.c.one
        zero = self.c.zero
        negone = self.c.negone

        if a.x.equal(zero) and a.y.equal(one):
            return self.set(self.c.s1)
        elif a.x.equal(zero) and a.y.equal(negone):
            return self.set(self.c.s2)
        elif a.x.equal(one) and a.y.equal(zero):
            return self.set(self.c.s3)
        elif a.x.equal(negone) and a.y.equal(zero):
            return self.set(self.c.s4)
        z = ModInt(self.c.p)
        self.x.set(a.y)
        self.y.set(a.x)
        self.z.set(z.mul(a.y, a.x))
        return self

    def to_ep(self, a):
        """
        Converts a point in inverted coordinates to its equivalent point in
        Edwards coordinates, with respect to the parent curve.
        """
        one = self.c.one
        zero = self.c.zero
        negone = self.c.negone
        ed = self.c.c.point()
        if a._special_pt(a):
            if a.equal(self.c.s1):
                ed.x.set(zero)
                ed.y.set(one)
            elif a.equal(self.c.s2):
                ed.x.set(zero)
                ed.y.set(negone)
            elif a.equal(self.c.s4):
                ed.x.set(negone)
                ed.y.set(zero)
            else:
                ed.x.set(one)
                ed.y.set(zero)
        else:
            x, y = ModInt(self.c.p), ModInt(self.c.p)
            x.div(a.z, a.x)
            y.div(a.z, a.y)
            ed.x.set(x)
            ed.y.set(y)
        return ed

Group.register(invEdwardsCurve)
Point.register(invEdwardsPoint)

Group.__subclasscheck__(invEdwardsCurve)
Point.__subclasscheck__(invEdwardsPoint)
Group.__instancecheck__(invEdwardsCurve)
Point.__instancecheck__(invEdwardsPoint)