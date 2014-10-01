import edwards
from modular import ModInt
from group import Group, Point

# TODO:
#   - conversions to weierstrauss form

class montEdwardsCurve(edwards.EdwardsCurve):
    """
    A Montgomery curve over a field K is defined by the equation:
    By^2 = x^3 + Ax^2 + X for A, B in K.
    K is a finite prime field of order p.

    Change of parameters maps:
    tE(a, d) --> M(A, B): (2(a + d)/(a - d), 4/(a - d))
    B = 4 * inverse(ed.a - ed.d, ed.p)
    A = 2 * (ed.a + ed.d) * inverse(ed.a - ed.d, ed.p)

    Curve25519 is the curve
    """

    def __init__(self, A, B, p, r, gx, gy, ed):
        self.name = "Montgomery"
        self.c = ed
        self.p = p
        self.r = r
        self.A = A
        self.B = B
        self.base = montEdwardsPoint(self, gx, gy)

        self.one = self.c.one
        self.zero = self.c.zero
        self.i = montEdwardsPoint(self.zero, self.one)

        if not self.base._on_curve():
            raise Exception("base is not on curve")

    def point(self):
        return montEdwardsPoint(self, ModInt(self.p), ModInt(self.p))

class montEdwardsPoint(edwards.EdwardsPoint):
    def __init__(self, curve, x=ModInt(), y=ModInt()):
        self.c = curve
        self.x = x
        self.y = y

    def _on_curve(self):
        """
        Tests that point (a) lies on the Montgomery curve with specified
        parameters: By^2 = x^3 + Ax^2 + x
        """
        x, y = self.x, self.y
        P = self.c.p
        l, r, yy, xx, t1 = ModInt(P), ModInt(P), ModInt(P), ModInt(P), ModInt(P)
        yy.mul(y, y)
        xx.mul(x, x)
        t1.mul(xx, self.c.A)
        l.mul(self.c.B, yy)
        r.exp(x, ModInt(P, 3)).add(r, t1).add(r, x)
        print("l", l.v, r.v, self)
        return l.equal(r) 

    def identity(self):
        self.set(self.c.i)
        return self

    def inverse(self, a):
        self.x.set(a.x)
        self.x.neg(a.y)
        return self

    def add(self, p, q):
        """
        Adds two points on the Montgomery curve together.
        Only certain points can be added together: p != +/- q
        So this addition law is not unified
        Computational cost:
        """
        print("p", p)
        x1, x2 = p
        y1, y2 = q

        if p == q or p == self.inverse(q):
            return self.double(p)

        if p == self.identity() and q == self.identity():
            return self.identity()
        elif p == self.identity():
            return q
        elif q == self.identity():
            return p

        delta = (y2 - y1) * inverse(x2 - x1, self.p)
        x3 = self.B * delta**2 - self.A - x1 - x2

        y3 = delta * (x1 - x3) - y1

        # firstnum = (2 * x1 + x2 + self.A) * (y2 - y1)
        # firstdenom = x2 - x1

        # secondnum = pow(y2 - y1, 3)
        # seconddenom = pow(x2 - x1, 3)
        # y3 = firstnum * inverse(firstdenom, self.p) - (self.B * secondnum) * inverse(seconddenom, self.p) - y1

        p = (x3 % self.p, y3 % self.p)
        print("\nincorrect answer", p)
        assert self.is_element(p)
        return p

    def double(self, p):
        """
        Doubles a point on the curve.
        Computational cost:
        """
        x1, y1 = p

        if p == self.identity():
            return self.identity()

        delta = (3 * x1**2 + 2 * self.A * x1 + 1) * inverse(2*self.B*y1, self.p)
        x3 = self.B * delta**2 - self.A - 2 *x1
        y3 = delta * (x1 - x3) - y1
        # num = (x1**2 - 1)**2
        # denom = (4*x1*(x1**2 + self.A * x1 + 1))
        # x3 = num * inverse(denom, self.p)

        # firstnum = (3 * x1 + self.A) * (3 * x1**2 + 2 * self.A*x1 +1)
        # firstdenom = 2 * self.B * y1

        # secondnum = self.B * pow(3 * x1**2 + 2 * self.A * x1 + 1, 3)
        # seconddenom = pow(2 * self.B * y1, 3)

        # y3 = firstnum * inverse(firstdenom, self.p) - secondnum * inverse(seconddenom, self.p) - y1

        p = (x3 % self.p, y3 % self.p)
        print("correct answer", p)
        assert self.is_element(p)
        return p


    def from_ep(self, a):
        return self

    def to_ep(self, a):
        return ed

