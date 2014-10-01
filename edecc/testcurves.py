from edwards import EdwardsCurve
import mont
from utils import string_to_long
from modular import ModInt

"""
These "safe curves" are obtained from
http://safecurves.cr.yp.to/ and
http://eprint.iacr.org/2013/647

Montgomery curves are defined by the following parameters:
    - p: order of the prime field
    - a: equation parameter
    - r: prime group order
    - (gx, gy): the base point, given by the x-coordinate and solving for the
    y-coordinate

(Twisted) Edwards curves are defined by the following parameters:
    - p: order of the prime field
    - d, a: equation parameters
    - r: prime group order
    - (gx, gy): the base point, given by the y coordinate and solving for the
    x-coordinate

Each curve includes a security level attribute, p-sec, as measured by the cost
of a Pollard p attack. For comparison, the NIST-certified P-256 curve has a
2^128 security level.

"""


# These curves haven't been tested against the arithmetic operations.

def m383():
    """ Returns the Montgomery Curve M-383, which has security level 2^189.8
    """
    p = pow(2, 383) - 187
    a = 2065150
    r = "2462625387274654950767440006258975862817483704404090416746934574\
        041288984234680883008327183083615266784870011007447"
    gx = 12
    gy = "47376234018917539976605463003759025768396171672577037256303897915244\
        63565757299203154901655432096558642117242906494"
    return mont.MontgomeryCurve(p, a, r, gx, gy)

def e222():
    """
    Returns an Edwards curve specially suited for crypto purposes: ED-222,
    discussed here:
    x^2+y^2 = 1+160102x^2y^2 ; modulo p = 2^222 - 117
    It has security level 2^109.8
    """
    p = pow(2, 222) - 117
    d = 160102
    s = 514330505974727608950741689135299334467771812183010135612886611057
    r = 1684996666696914987166688442938726735569737456760058294185521417407
    gx = 2705691079882681090389589001251962954446177367541711474502428610129
    gy = 28
    return edecc.EdwardsCurve(p, d, 1, r, gx, gy)

def numsp256t1():
    p = string_to_long("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\
        FFFFFFFFFFFFF43")
    a = string_to_long("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\
        FFFFFFFFFFFF42")
    d = string_to_long("3BEE")
    r = string_to_long("3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBE6AA55AD0A6BC64E5B\
        84E6F1122B4AD")
    gx = 13
    gy = 56558017671127616723941111350859963393592894419123286895630403139519716830650
    return edecc.EdwardsCurve(p, d, 1, r, gx, gy)

def e382():
    """ Returns the Edwards Curve E-382, which has security level 2^189.8
    """
    p = 2^382 - 105
    d = -67254
    r = 2462625387274654950767440006258975862817483704404090416745738034557663054564649171262659326683244604346084081047321
    gx = 3914921414754292646847594472454013487047137431784830634731377862923477302047857640522480241298429278603678181725699
    gy = 13
    return edecc.EdwardsCurve(p, d, 1, r, gx, gy)

def curve1174():
    """ Returns the Edwards Curve curve1174, which has security level 2^124.3
    """
    p = pow(2, 251) - 9
    d = -1174
    r = "9046256971665327767466483203803742800923390352794954740234892617\
        73642975601"
    gx = "1582619097725911541954547006453739763381091388846394833492296309\
        729998839514"
    gy = "30375380136041545047641157286514376465195135343052234227548270556\
        89195992590"
    return edecc.EdwardsCurve(p, d, 1, r, gx, gy)

