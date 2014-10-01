#!/usr/bin/env python
from Crypto.Hash import SHA256
from Crypto.Util.number import bytes_to_long, GCD

class ElGamal:
    def _hash(self, data, bits):
        limit = bits // 8
        h = SHA256.new()
        h.update(data)
        return bytes_to_long(h.digest()[0:limit])

    def encrypt(self, element, data):
        y = self.secret()
        s = self.point()
        s.multiply(element, y)

        c1 = self.point()
        c1.multiply(element.generator(), y)

        de = self.point().encode(data)
        c2 = self.point().add(de, s)

        return (c1, c2)

    def decrypt(self, secret, encrypted):
        c1, c2 = encrypted
        s = self.point()
        s.multiply(c1, secret)

        de = self.point()
        de.sub(c2, s)
        return de.decode(de)

class PublicKey:
    def __init__(self, group, element):
        self.group = group
        self.element = element

    def exchange(self, private):
        pubkey = self.group.point()
        if isinstance(self.element.string(), tuple):
            return pubkey.multiply(self.element, private.secret).x.string()
        return pubkey.multiply(self.element, private.secret)

    def verify(self, data, sign):
        return self.group.verify(self.element, data, sign)

    def encrypt(self, data):
        return self.group.encrypt(self.element, data)

class PrivateKey(PublicKey):
    def __init__(self, group, secret = None):
        if secret == None:
            secret = group.secret()
        element = group.point()
        element.random_element()
        super(PrivateKey, self).__init__(group, element)
        self.secret = secret

    def exchange(self, public):
        privatekey = self.group.point()
        if isinstance(self.element.string(), tuple):
            return privatekey.multiply(public.element, self.secret).x.string()
        return privatekey.multiply(public.element, self.secret)

    def sign(self, data):
        return self.group.sign(self.secret, data)

    def decrypt(self, encrypted):
        return self.group.decrypt(self.secret, encrypted)

    def public_key(self):
        return PublicKey(self.group, self.element)
