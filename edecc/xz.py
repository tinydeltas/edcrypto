import edwards
from modular import ModInt
from group import Group, Point

class xzEdwardsCurve(edwards.EdwardsCurve):
    """
    Alternate coordinate system for Montgomery curves: XZ (projective) coordinates
    M(x, y) is represented as x = X/Z
     By^2 = x^3 + Ax^2 + X

    The xz-coordinate system allows for fast x-coordinate computation for
    Montgomery curves.

    Originally introduced in 1987 by Montgomery in this paper:
    http://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866113-7/S0025-5718-1987-0866113-7.pdf

    """

    def __init__(self, mont):
        self.name = "Montgomery xz"
        self.c = mont
        self.a = mont.A
        #self.d = mont.d
        self.p = mont.p
        self.r = mont.r

        self.zero = mont.zero
        self.one = mont.one

        self.base = self.point().from_mp(mont.base)
        self.i = xzEdwardsPoint(self, self.zero, self.zero) # uncertain
        if not self.base._on_curve():
            raise Exception("Incorrect base point")

    def point(self):
        return xzEdwardsPoint(self, ModInt(self.p), ModInt(self.p))

class xzEdwardsPoint(edwards.EdwardsPoint):
    def __init__(self, curve, x=ModInt(), z=ModInt()):
        self.c = curve
        self.x = x
        self.z = z

    def string(self):
        return (self.x.v, self.z.v)

    def _on_curve(self):
        """
        Tests that point (a) lies on the Montgomery curve with specified
        parameters: By^2 = x^3 + Ax^2 + x
        """
        print("P", self.c.p)
        newpt = self.to_mp(self)
        return super(newpt)._on_curve()

    def identity(self):
        self.set(self.c.i)
        return self

    def inverse(self, a):
        self.x.set(a.x)
        self.z.neg(a.z)
        return self

    def equal(self, p):
        return self.x.equal(p.x) and self.z.equal(p.z)

    def set(self, p):
        self.x.set(p.x)
        self.z.set(p.z)
        return self

    def add(self, p, q, m):
        """
        Addition in XZ coordinates works as follows:

        """
        xp, zp = p
        xq, zq = q
        xm, zm = m

        x = 4 * (xq * xp - zq * zp)**2 * zm
        y = 4 * (xq * zp - zq * xp)**2 * xm
        return (x % self.p, y % self.p)

    def double(self, p):
        """
        Doubling in XZ coordinates.
        Computational cost: 2M + 2S = 4M (3M if z == 1), "dbl-1987-m-3"
        """
        X, Z = p
        assert X != 0 or Z != 0
        if Z == 0:
            return self.identity()

        if Z == 1:
            XX = X**2
            x_r = (XX - 1)**2
            y_r = 4 * X * (XX + self.A * X + 1)
        else:
            A = X + Z
            AA = A**2
            B = X - Z
            BB = B**2
            C = AA - BB
            x_r = AA * BB
            z_r = C * (BB + self.A * 25 * C)
        return (x_r % self.p, z_r % self.p)

    def multiply(self, n, P):
        """
        Source: http://cr.yp.to/ecdh/curve25519-20051115.pdf
        """
        if n == 1:
            return P
        elif n == 2:
            return self.double(P)

        if n % 2 == 0:
            return 2 * self.multiply(n//2, P)
        else:
            first = self.add(self.multiply(n//2, P), P)
            return first + self.multiply(n//2, P)

    def yrecover(self, x):
        """
        Sanity checker.
        y^2 = (x^3 +Ax^2 + x) / B
        """
        print("current group", self.c)
        yy = ModInt(p=self.c.p)
        yy.exp(x, ModInt(v=3)).add(yy.mul(self.c.A, yy.mul(x, x)), x).div(yy, self.c.B)
        if yy.jacobi(yy) == 1:
            return yy.sqrt(yy)

    def from_mp(self, a):
        self.x.set(a.x)
        self.z.set(self.c.one)
        return self

    def to_mp(self, a):
        return edwards.EdwardsPoint()

Group.register(xzEdwardsCurve)
Point.register(xzEdwardsPoint)

Group.__subclasscheck__(xzEdwardsCurve)
Point.__subclasscheck__(xzEdwardsPoint)
Group.__instancecheck__(xzEdwardsCurve)
Point.__instancecheck__(xzEdwardsPoint)