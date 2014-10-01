from mont import montEdwardsCurve
from modular import ModInt
from ed25519 import ed25519

def curve25519():
    """
    Returns an instance of Curve25519, a Montgomery curve described
    by Daniel J. Bernstein in his paper Curve25519.

    Curve25519 is of the form
    y^2 = x^3 + 486662x^2 + x, with associated p and base.
    Its security level is 2^125.8

    Some of its advantages include:
        - having the prime p be congruent to 5 mod 8, making the square root
        calculation easier
        - p is of the form p = 2^m - d, where d < ceil(log(p)) == m
    """
    prime = pow(2, 255) - 19
    p = ModInt(prime, prime)
    A = ModInt(p, 486662)
    B = ModInt(p, 1)
    r = ModInt(p, 7237005577332262213973186563042994240857116359379907606001950938285454250989)
    gx = ModInt(p, 9)
    gy = ModInt(p, 14781619447589544791020593568409986887264606134616475288964881837755586237401)
    return montEdwardsCurve(A, B, p, r, gx, gy, ed25519())