import edwards
import proj
import inv
import mont
import extended
import xz
from modular import ModInt

def pp_to_ep(curve, a):
        """
        Converts a, a point in projective Edwards coordinates, to a point in
        standard Edwards coordinates.
        """
        return ((a[0] * inverse(a[2], self.p)) % self.p,
            (a[1] * inverse(a[2], self.p)) % self.p)

def ep_to_pp(curve, a):
    return proj.projEdwardsPoint(a.c, a.x, a.y, a.c.one)

def pp_to_ip(curve, a):
    """
    Converts a, a point in projective Edwards coordinates, to a point in
    inverted Edwards coordinates.
    """
    return ((a[1]*a[2]) % self.p, (a[0]*a[2]) % self.p, (a[0]*a[1]) % self.p)

def pp_to_exp(curve, a):
    """
    Converts a, a point in projective Edwards coordinates, to a point in
    extended Twisted Edwards coordinates.
    """
    X, Y, Z = a
    return ((X*Z) % self.p, (Y*Z) % self.p, Z, Z**2 % self.p)

# Inverted

def ep_to_ip(c, a):
    # I'm sorry this is really disgusting
    negone = ModInt(p=c.p, v=-1)
    if (a.x.v, a.y.v) == (0, 1):
        return invEdwardsPoint(c, c.one, c.zero, c.zero)
    elif (a.x.v, a.y.v) == (0, -1):
        return invEdwardsPoint(c, negone, c.zero, c.zero)
    elif (a.x.v, a.y.v) == (1, 0):
        return invEdwardsPoint(c, c.zero, negone, c.zero)
    elif (a.x.v, a.y.v) == (-1, 0):
        return invEdwardsPoint(c, c.zero, c.one, c.zero)
    else:
        z = ModInt(p=c.p)
        return inv.invEdwardsPoint(a.c, a.y, a.x, z.mul(a.y, a.x))

def ip_to_ep(curve, a):
        """
        Converts a, a point in inverted Edwards coordinates, to a point in
        standard Edwards coordinates.
        """
        if a == (1, 0, 0):
            return (0, 1)
        elif a == (-1, 0, 0):
            return (0, -1)
        elif a == (0, 1, 0):
            return (-1, 0)
        elif a == (0, -1, 0):
            return (1, 0)
        return ((a[2] * inverse(a[0], self.p)) % self.p,
                (a[2] * inverse(a[1], self.p)) % self.p)

#Extended
#
def exp_to_ep(a):
    """
    Converts a, a point in extended Edwards coordinates, to a point in
    standard Edwards coordinates.
    Alternative method: ext -> projective -> regular
    """
    X = a.x; Y = a.y; T = a.t; Z = a.z
    x = ModInt(p=a.c.p); y = ModInt(p=a.c.p)
    x.div(X, Z)
    y.div(Y, Z)

    t1 = ModInt(p=a.c.p); t2 = ModInt(p=a.c.p)
    assert t1.mul(x, y).equal(t2.div(T, Z))
    return edwards.EdwardsPoint(a.c, x, y)


def ext_to_pp(self, a):
    """
    Converts a, a point in extended Edwards coordinates, to a point in
    projective Edwards coordinates.
    """
    return (a[0], a[1], a[3])


def ext_to_exc(c):
    """
    Uses map (x, y) -> (x/sqrt(-a), y) to define more compact extended
    curve.
    """
    pass

def to_ec(c):
    return c.c

# Montgomery

def mp_to_ep(a):
    """
    Converts a, a point in Montgomery coordinates, to
    a point on a (twisted) Edwards curve.
    """
    u, v = a
    x = u * inverse(v, self.p)
    y = (u - 1) * inverse((u + 1), self.p)
    p = (x % self.p, y % self.p)
    print(p)
    return p

def mp_to_xz(a):
    """
    Converts a, a point in Montgomery coordinates, to
    XZ (Montgomery) coordinates
    """
    print("p?", a.c.p)
    return xz.xzEdwardsPoint(a.c, a.x, a.c.one)

def mc_to_tc(c):
    """
    Returns a twisted Edwards curve from Montgomery curve parameters.
    """
    d = (self.A - 2) * inverse(self.B, self.p)
    a = (self.A + 2) * inverse(self.B, self.p)
    return EdwardsCurve(self.p, d, a, self.r, self.base[0], self.base[1])

def mc_to_wc(c):
    """
    Converts curve c (Montgomery) to an equivalent curve in Weierstrauss
    form.
    """
    a = (3 - self.A**2) * inverse(3 * self.B**2, self.p)
    b = (2 * self.A**3 - 9 * self.A) * inverse(27 * self.B**3, self.p)
    k = ""
    gx = ""
    gy = ""
    return EllipticCurve(a % self.p, b % self.p, self.p, self.r, k, gx, gy)

def xz_to_mp(a):
    x = ModInt(p=a.c.p)
    x.div(a.x, a.z)
    y = ModInt(p=a.c.p)
    y = a.yrecover(x)
    return mont.montEdwardsPoint(a.c, x, y)

def xz_to_mp(a):
    """
    Converts a, a point in XZ - Montgomery coordinates, to
    a point on the equivalent Montgomery curve
    """
    x = ModInt(p=a.c.p)
    x.div(a.x, a.z)
    y = ModInt(p=a.c.p)
    return mont.montEdwardsPoint(a.c, x, a.yrecover(x))

def to_tp(self, a):
        """
        Converts a, a point in standard Edwards coordinates, to
        twisted Edwards coordinates.
        (x, y) --> (x, y/a)
        """
        self.x = a.x
        self.y = ModInt().div(a.y, self.c.a)
        return self

def from_tp_to_ep(self, a):
    """
    Converts a, a point in twisted Edwards coordinates, to
    standard Edwards coordinates.
    (x, y) --> (x/sqrt(a), y)
    """
    x = ModInt(); sq = ModInt()
    self.x = x.div(a.x, sq.sqrt(self.c.a))
    self.y = a.y
    return self

def from_ep(self, a):
    return a

def to_ep(self, a):
    return a

def ep_to_mp(self, a):
    """
    Converts a, a point in twisted Edwards coordinates, to
    Montgomery coordinates.
    """
    x = a.x
    y = a.y
    u = (1 + y) * inverse(1 - y, self.p)
    v = -(1 + y) * (inverse((1 - y) * x, self.p))
    print("v", self.p - v % self.p)
    return (u % self.p, v % self.p)
