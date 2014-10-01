"""
Public-key cryptography abstract-group interface.
"""

from abc import ABCMeta, abstractmethod

class Group(object, metaclass=ABCMeta):
    """
    Abstract cryptographic group suitable as the underlying secure
    mathematical structure to schemes such as Diffie-Hellman key exchange,
    ElGamal, etc.
    """

    @abstractmethod
    def string(self):
        """
        Returns name of the group.
        """
        pass

    @abstractmethod
    def order(self):
        """
        Returns order of the group.
        """
        pass

    @abstractmethod
    def point(self):
        """
        Create new element within the group.
        """
        pass

    @abstractmethod
    def secret(self):
        """
        Create new secret (scalar) within the group.
        """
        pass

class Secret(object, metaclass=ABCMeta):
    """"
    A secret represents a scalar, used in encrypting an element of some
    cryptographically secure group; e.g., exponent in DSA finite fields
    (relying on DLP), scalar multiplier in elliptic curve groups.

    This class contains modular arithmetic methods as well as limited
    interfacing with the group.

    Attributes:
        - v: the value represented by the secret.
        - p: as in Z_p, the modulo.
    """

    @abstractmethod
    def zero(self):
        """
        The zero element (additive identity).
        """
        pass

    @abstractmethod
    def one(self):
        """
        The one element (multiplicative identity).
        """
        pass

    @abstractmethod
    def add(self, a, b):
        """
        Modular addition.
        """
        pass

    @abstractmethod
    def sub(self, a, b):
        """
        Modular subtraction.
        """
        pass

    @abstractmethod
    def neg(self, a):
        """
        Modular negation.
        """
        pass

    @abstractmethod
    def mul(self, a, b):
        """
        Modular product.
        """
        pass

    @abstractmethod
    def div(self, a):
        """
        Modular division.
        """
        pass

    @abstractmethod
    def inv(self, a):
        """
        Modular inverse.
        """
        pass

    @abstractmethod
    def equal(self, a):
        pass

    @abstractmethod
    def set(self, a):
        """
        Set the properties of the current ModInt object(self) to those of a.
        """
        pass

    @abstractmethod
    def random_secret(self):
        """
        PIcks a random secret in the range [1, modulo p).
        """
        pass

class Point(object, metaclass=ABCMeta):
    """
    Represents an element of a public-key cryptographic group, e.g. ,
    number modulo prime for Schnorr groups and finite fields,
    (x, y) point on an elliptic key.
    """

    @abstractmethod
    def _on_curve(self):
        """
        Checks if the point is on the curve defined by its parent group.
        """
        pass

    @abstractmethod
    def generator(self):
        """
        Sets the point to the base point, or generator, of its parent group
        and returns it.
        """
        pass

    @abstractmethod
    def identity(self):
        """
        Sets the point to the additive identity of its parent group
        and returns it.
        """
        pass

    @abstractmethod
    def neg(self, a):
        """
        Negates a according to the negation law of its parent group,
        sets the point to -a, and returns it.
        """
        pass

    @abstractmethod
    def equal(self, a):
        """
        Checks if the point equals the point a.
        """
        pass

    @abstractmethod
    def set(self, a):
        """
        Sets the point to a by setting each of its individual components.
        """
        pass

    @abstractmethod
    def add(self, a, b):
        """
        Adds the points a and b according to the addition operation of its
        parent group, sets the point to the sum, and returns it.
        """
        pass

    @abstractmethod
    def sub(self, a, b):
        """
        Adds the points a and neg(b) according to the addition operation of its
        parent group, sets the point to the sum, and returns it.
        """
        pass

    @abstractmethod
    def double(self, a):
        """
        Doubles the point a, sets the point to the sum, and returns it.
        """
        pass

    @abstractmethod
    def multiply(self, P, n):
        """
        Multiplies the point P by the scalar n, sets the point to the product,
        and returns it.
        """
        pass

    @abstractmethod
    def encode(self, data):
        """
        Encodes a limited amount of data in a randomly chosen point.
        """
        pass

    @abstractmethod
    def decode(self, pt):
        """
        Decodes data contained in the point.
        """
        pass

    @abstractmethod
    def random_element(self, p):
        """
        Picks a random element from the group.
        """
        pass

