import edwards
from utils import string_to_long
from modular import ModInt

def ed25519():
    """ Edwards Curve version of curve25519.
    Base points obtained from
    http://tools.ietf.org/html/draft-ladd-safecurves-04
    """
    prime = pow(2, 255) - 19
    p = ModInt(prime, prime)
    d = ModInt(p)
    d.div(ModInt(p, -121665), ModInt(p, 121666))
    a = ModInt(p, -1)
    r = ModInt(p, (pow(2, 252) + 27742317777372353535851937790883648493 % p.v))
    gx = ModInt(p, string_to_long("216936 d3cd6e 53fec0 a4e231 fdd6dc 5c692c c76095 25a7b2 c9562d 608f25 d51a"))
    gy = ModInt(p)
    gy.div(ModInt(p, 4), ModInt(p, 5))
    return edwards.EdwardsCurve("Twisted Edwards w/ a = -1", p, d, a, r, gx, gy)