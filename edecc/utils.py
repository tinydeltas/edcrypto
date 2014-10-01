"""
Misc helper functions.
"""

import binascii
from Crypto.Util.number import bytes_to_long, long_to_bytes

def b2l(s):
    return bytes_to_long(s)

def l2b(s):
    return long_to_bytes(s)

def string_to_long(s):
    s = bytes("".join(s.split()), "UTF-8")
    s = binascii.a2b_hex(s)
    return bytes_to_long(s)
